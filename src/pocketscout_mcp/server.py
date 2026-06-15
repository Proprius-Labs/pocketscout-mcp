"""PocketScout MCP Server

Scout the binding landscape before you design the binder.

A FastMCP server for drug-target binding-site triage: it aggregates
structural, chemical, conservation, and literature data into a briefing
on a protein's druggable pockets. Especially useful as the reconnaissance
step before computational binder design — the gap between "I have a target"
and "I'm running RFdiffusion."

Tools are ordered to reflect the natural scientific workflow:
1. CharacterizeTarget — biological context + AlphaFold confidence (per-residue pLDDT)
2. GetRelatedStructures — structural coverage for the target
3. GetBindingSites — map known pockets from co-crystals (with contact types)
4. GetLigandHistory — competitive landscape from ChEMBL
5. CheckConservation — human vs. mouse/rat/cynomolgus translatability
6. SearchTargetLiterature — recent structural/design insights
7. CheckKnownVariants — disease/resistance variants at binding residues
8. ConsolidateBindingSites — union of pockets across all structures

Two prompts orchestrate the tools: target_briefing (quick triage) and
binding_site_assessment (in-depth ranked recommendation).
"""

from __future__ import annotations

from contextlib import asynccontextmanager

from fastmcp import FastMCP
from mcp.types import Icon
from starlette.requests import Request
from starlette.responses import Response

from .clients.uniprot import UniProtClient, parse_target_profile, parse_known_variants
from .clients.pdb import PDBClient, parse_structure_metadata, ARTIFACT_LIGANDS
from .clients.alphafold import AlphaFoldClient, analyze_confidence
from .clients.chembl import ChEMBLClient, analyze_bioactivities
from .clients.pubmed import PubMedClient
from .clients.base import APIError
from .models import (
    TargetProfile,
    TopologyRegion,
    ConfidenceRegion,
    ConstructCoverage,
    BindingSite,
    BindingSiteResidue,
    BindingSiteMap,
    RelatedStructure,
    RelatedStructuresResult,
    LigandHistory,
    ResidueConservation,
    SpeciesConservation,
    ConservationResult,
    KnownVariant,
    VariantCheckResult,
    PaperResult,
    LiteratureResult,
    ConsolidatedPocket,
    ConsolidationResult,
)

# ---------------------------------------------------------------------------
# Server + clients
# ---------------------------------------------------------------------------

ICON_URL = "https://raw.githubusercontent.com/Proprius-Labs/pocketscout-mcp/main/assets/icon.svg"
ICONS = [Icon(src=ICON_URL, mimeType="image/svg+xml")]

ICON_SVG = b"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round" role="img" aria-label="PocketScout">
  <style>
    svg { color: #1a1a1a; }
    @media (prefers-color-scheme: dark) {
      svg { color: #f5f5f5; }
    }
  </style>
  <line x1="12" y1="2.5" x2="12" y2="5"/>
  <line x1="12" y1="19" x2="12" y2="23"/>
  <line x1="1" y1="12" x2="6" y2="12"/>
  <line x1="18" y1="12" x2="23" y2="12"/>
  <polygon points="12,0.5 10.5,2.5 13.5,2.5" fill="currentColor" stroke="none"/>
  <polygon points="12,5 18,8.5 18,15.5 12,19 6,15.5 6,8.5"/>
  <circle cx="18" cy="8.5" r="0.7" fill="currentColor" stroke="none"/>
  <circle cx="12" cy="19" r="0.7" fill="currentColor" stroke="none"/>
  <circle cx="6" cy="15.5" r="0.7" fill="currentColor" stroke="none"/>
  <circle cx="6" cy="8.5" r="0.7" fill="currentColor" stroke="none"/>
  <circle cx="12" cy="5" r="1.2" fill="#0D9488" stroke="none"/>
  <circle cx="18" cy="15.5" r="1.2" fill="#0D9488" stroke="none"/>
  <circle cx="12" cy="12" r="2" stroke="#D97706" stroke-width="0.6" fill="none"/>
  <circle cx="12" cy="12" r="0.95" fill="#D97706" stroke="none"/>
</svg>"""

uniprot = UniProtClient()
pdb = PDBClient()
alphafold = AlphaFoldClient()
chembl = ChEMBLClient()
pubmed = PubMedClient()


@asynccontextmanager
async def _lifespan(app):
    try:
        yield
    finally:
        for c in (uniprot, pdb, alphafold, chembl, pubmed):
            await c.close()


mcp = FastMCP(
    "PocketScout",
    instructions=(
        "Fast binding-site intelligence for drug-target triage. Pulls together "
        "structural, chemical, conservation, and literature data into a briefing "
        "on a protein's druggable pockets. Use the target_briefing prompt for a "
        "quick assessment, or binding_site_assessment for an in-depth, "
        "design-focused workup."
    ),
    icons=ICONS,
    lifespan=_lifespan,
)


@mcp.custom_route("/favicon.ico", methods=["GET"], include_in_schema=False)
@mcp.custom_route("/favicon.svg", methods=["GET"], include_in_schema=False)
async def favicon(request: Request) -> Response:
    return Response(
        content=ICON_SVG,
        media_type="image/svg+xml",
        headers={"Cache-Control": "public, max-age=86400"},
    )


# ---------------------------------------------------------------------------
# Helper: resolve PDB ID → UniProt accession
# ---------------------------------------------------------------------------

async def _resolve_identifiers(pdb_id: str | None = None, uniprot_id: str | None = None) -> tuple[str | None, str | None]:
    """Resolve PDB ID and UniProt accession from either input.

    Many tools need both identifiers. This helper does the cross-reference
    so individual tools don't have to duplicate the logic.
    """
    if pdb_id and not uniprot_id:
        mappings = await pdb.get_uniprot_mapping(pdb_id)
        uniprot_id = mappings[0] if mappings else None
    return pdb_id, uniprot_id


# ---------------------------------------------------------------------------
# Tool 1: characterize_target
# ---------------------------------------------------------------------------

READ_ONLY_ANNOTATIONS = {
    "readOnlyHint": True,
    "destructiveHint": False,
    "idempotentHint": True,
    "openWorldHint": True,
}


@mcp.tool(
    name="CharacterizeTarget",
    title="Characterize Target",
    tags={"target-profile", "structure"},
    annotations=READ_ONLY_ANNOTATIONS,
)
async def characterize_target(
    pdb_id: str | None = None,
    uniprot_id: str | None = None,
) -> dict:
    """Establish biological context for a drug target protein.

    Retrieves protein function, family classification, subcellular location,
    disease associations, and AlphaFold structure confidence. This should
    be your FIRST call — all downstream analysis depends on this context.

    IMPORTANT: AlphaFold confidence flags regions where predicted structure
    is unreliable. Low-confidence regions (<70 pLDDT) may have incorrect
    pocket predictions — always cross-reference with experimental structures.

    Provide EITHER a PDB ID (e.g. '1M17') OR a UniProt accession (e.g. 'P00533').
    If a PDB ID is given, the UniProt mapping is resolved automatically.
    """
    if not pdb_id and not uniprot_id:
        return {"error": "Provide either pdb_id or uniprot_id"}

    pdb_id, uniprot_id = await _resolve_identifiers(pdb_id, uniprot_id)

    if not uniprot_id:
        return {"error": f"Could not map PDB {pdb_id} to a UniProt accession. The structure may contain only non-standard chains."}

    # Fetch UniProt data
    try:
        entry = await uniprot.get_entry(uniprot_id)
    except APIError as e:
        return {"error": f"Failed to fetch UniProt entry: {e}"}

    profile_data = parse_target_profile(entry)

    # Fetch AlphaFold confidence
    try:
        af_prediction = await alphafold.get_prediction(uniprot_id)
        try:
            per_residue = await alphafold.get_per_residue_plddt(uniprot_id)
        except APIError:
            per_residue = []
        af_data = analyze_confidence(af_prediction, profile_data["sequence_length"], per_residue or None)
    except APIError:
        af_data = {
            "overall_confidence": None,
            "regions": [],
            "warnings": ["Could not reach AlphaFold DB."],
        }

    # Build confidence regions
    confidence_regions = [
        ConfidenceRegion(**r) for r in af_data.get("regions", [])
    ]

    # Build topology regions
    topology_regions = [
        TopologyRegion(**t) for t in profile_data.pop("topology", [])
    ]

    profile = TargetProfile(
        **profile_data,
        topology=topology_regions,
        alphafold_confidence=confidence_regions,
        alphafold_overall_confidence=af_data.get("overall_confidence"),
        low_confidence_warnings=af_data.get("warnings", []),
    )

    return profile.model_dump()


# ---------------------------------------------------------------------------
# Tool 2: get_related_structures
# ---------------------------------------------------------------------------

@mcp.tool(
    name="GetRelatedStructures",
    title="Get Related Structures",
    tags={"structure"},
    annotations=READ_ONLY_ANNOTATIONS,
)
async def get_related_structures(
    pdb_id: str | None = None,
    uniprot_id: str | None = None,
    limit: int = 20,
) -> dict:
    """Find all PDB structures for a target protein.

    Returns all available experimental structures, sorted by resolution.
    Use this to understand structural coverage: how many structures exist,
    what ligands have been co-crystallized, what conformational states
    are captured, and what the best-quality structure is.

    A target with many high-resolution co-crystal structures has a rich
    binding site landscape to analyze. A target with only 1-2 structures
    (or only apo/unliganded structures) has less structural evidence.

    Call this AFTER CharacterizeTarget and BEFORE GetBindingSites to
    identify which structures to analyze for pockets.
    """
    if not pdb_id and not uniprot_id:
        return {"error": "Provide either pdb_id or uniprot_id"}

    pdb_id, uniprot_id = await _resolve_identifiers(pdb_id, uniprot_id)

    if not uniprot_id:
        return {"error": f"Could not resolve UniProt ID from PDB {pdb_id}"}

    # Search for all structures
    try:
        pdb_ids = await pdb.search_by_uniprot(uniprot_id, limit=limit)
    except APIError as e:
        return {"error": f"PDB search failed: {e}"}

    # Fetch UniProt profile once — provides gene name and organism for all structures
    organism = "unknown"
    gene_name = ""
    try:
        profile = parse_target_profile(await uniprot.get_entry(uniprot_id))
        gene_name = profile["gene_name"]
        organism = profile["organism"]
    except APIError:
        pass

    if not pdb_ids:
        return {
            "uniprot_id": uniprot_id,
            "gene_name": gene_name,
            "total_structures": 0,
            "structures": [],
            "unique_ligands": [],
            "has_apo_structure": False,
            "best_resolution": None,
            "interpretation": "No experimental structures found in PDB. Consider using AlphaFold predicted structure, but note that predicted binding sites have lower reliability than those from experimental structures.",
        }

    # Fetch metadata for each structure in parallel with bounded concurrency
    import asyncio
    sem = asyncio.Semaphore(8)

    async def _fetch(pid: str) -> RelatedStructure | None:
        async with sem:
            try:
                entry = await pdb.get_entry(pid)
                meta = parse_structure_metadata(entry)
                np_entities = await pdb.get_nonpolymer_entities(pid)
            except APIError:
                return None
        ligand_ids = []
        for np in np_entities:
            ids = np.get("rcsb_nonpolymer_entity_container_identifiers", {})
            comp_id = ids.get("nonpolymer_comp_id", "") or ids.get("comp_id", "")
            if not comp_id:
                pdbx = np.get("pdbx_entity_nonpoly", {})
                if isinstance(pdbx, list) and pdbx:
                    comp_id = pdbx[0].get("comp_id", "")
                elif isinstance(pdbx, dict):
                    comp_id = pdbx.get("comp_id", "")
            if comp_id and comp_id not in ARTIFACT_LIGANDS:
                ligand_ids.append(comp_id)
        return RelatedStructure(
            pdb_id=pid,
            title=meta["title"],
            resolution=meta["resolution"],
            method=meta["method"],
            organism=organism,
            ligand_ids=ligand_ids,
            release_date=meta["release_date"],
        )

    fetched = await asyncio.gather(*[_fetch(pid) for pid in pdb_ids[:limit]])
    structures = [s for s in fetched if s is not None]
    all_ligands = set()
    has_apo = False
    for s in structures:
        all_ligands.update(s.ligand_ids)
        if not s.ligand_ids:
            has_apo = True

    # Sort by resolution
    structures.sort(key=lambda s: s.resolution if s.resolution else 99.0)
    best_res = structures[0].resolution if structures else None

    # Interpretation
    interp_parts = [f"{len(structures)} experimental structures found for {gene_name or uniprot_id}."]
    if best_res:
        interp_parts.append(f"Best resolution: {best_res:.1f} Å.")
    interp_parts.append(f"{len(all_ligands)} unique ligands co-crystallized.")
    if has_apo:
        interp_parts.append("Apo (unliganded) structure(s) available for conformational comparison.")
    if len(structures) > 10:
        interp_parts.append("Rich structural coverage — multiple conformations and ligand-bound states available for comprehensive binding site analysis.")
    elif len(structures) < 3:
        interp_parts.append("Limited structural data. Binding site analysis will be constrained to available structures.")

    result = RelatedStructuresResult(
        uniprot_id=uniprot_id,
        gene_name=gene_name,
        total_structures=len(structures),
        structures=structures,
        unique_ligands=sorted(all_ligands),
        has_apo_structure=has_apo,
        best_resolution=best_res,
        interpretation=" ".join(interp_parts),
    )
    return result.model_dump()


# ---------------------------------------------------------------------------
# Tool 3: get_binding_sites
# ---------------------------------------------------------------------------

@mcp.tool(
    name="GetBindingSites",
    title="Get Binding Sites",
    tags={"structure", "binding-site"},
    annotations=READ_ONLY_ANNOTATIONS,
)
async def get_binding_sites(pdb_id: str) -> dict:
    """Map all known binding sites in a protein structure from co-crystallized ligands.

    Identifies binding pockets by analyzing non-polymer entities (ligands,
    cofactors) in the structure, filtering out crystallization artifacts
    (glycerol, PEG, sulfate, etc.), and classifying each site by type.

    Each site includes druggability assessment and modality recommendations:
    - Small deep pockets (< 20 contact residues) favor small molecules
    - Large flat interfaces (> 30 residues) favor biologics or de novo protein binders
    - Allosteric sites may offer selectivity advantages over orthosteric sites

    Call this on specific PDB IDs identified by GetRelatedStructures.
    For comprehensive analysis, call on multiple structures with different
    co-crystallized ligands to build a complete pocket map.
    """
    try:
        site_data = await pdb.get_binding_sites(pdb_id)
    except APIError as e:
        return {"error": f"Failed to fetch binding site data: {e}"}

    ligands = site_data.get("ligands", [])
    artifacts = site_data.get("artifacts_filtered", [])

    if not ligands:
        return BindingSiteMap(
            pdb_id=pdb_id.upper(),
            num_sites=0,
            sites=[],
            artifact_ligands_filtered=artifacts,
            interpretation=f"No non-artifact ligands found in {pdb_id}. This is an apo structure or contains only crystallization artifacts ({', '.join(artifacts) if artifacts else 'none detected'}). Check other structures for this target via GetRelatedStructures.",
        ).model_dump()

    # Build binding sites from ligand data
    sites = []
    for i, lig in enumerate(ligands):
        comp_id = lig.get("comp_id", "unknown")
        name = lig.get("name", "")

        # Classify site type based on ligand properties
        site_type = _classify_site_type(comp_id, name)

        # Build residue contacts from coordinate analysis
        contacts = lig.get("contacts", [])
        residues = [
            BindingSiteResidue(
                chain=c["chain"],
                residue_name=c["residue_name"],
                residue_number=c["residue_number"],
                contact_type=c.get("contact_type", "unknown"),
            )
            for c in contacts
        ]
        residue_positions = sorted({c["residue_number"] for c in contacts if c["residue_number"] > 0})

        site = BindingSite(
            site_id=f"site_{i+1}_{comp_id}",
            residues=residues,
            residue_positions=residue_positions,
            ligand_id=comp_id,
            ligand_name=name,
            site_type=site_type,
            evidence="co-crystal",
            num_residues=len(residues),
            druggability_notes=_assess_druggability(site_type, comp_id, len(residues)),
        )
        sites.append(site)

    # Reclassify sites: if a structure has both a cofactor/substrate site and
    # another site with non-overlapping residues, the second site is likely allosteric
    if len(sites) > 1:
        cofactor_positions = set()
        for s in sites:
            if s.site_type == "cofactor":
                cofactor_positions.update(s.residue_positions)
        if cofactor_positions:
            for s in sites:
                if s.site_type == "orthosteric" and s.residue_positions:
                    overlap = cofactor_positions.intersection(s.residue_positions)
                    overlap_frac = len(overlap) / len(s.residue_positions) if s.residue_positions else 1.0
                    if overlap_frac < 0.3:
                        s.site_type = "allosteric"
                        s.druggability_notes = _assess_druggability("allosteric", s.ligand_id or "", s.num_residues)

    # ---- Construct coverage check ----
    coverage = await _compute_construct_coverage(pdb_id)

    # ---- Build interpretation ----
    interp_parts = []

    # Lead with accessibility warning if present — this is the most important finding
    if coverage and coverage.accessibility_warning:
        interp_parts.append(f"⚠ ACCESSIBILITY: {coverage.accessibility_warning}")

    interp_parts.append(f"{len(sites)} binding site(s) identified from co-crystallized ligands in {pdb_id}.")

    if coverage:
        interp_parts.append(
            f"Structure covers residues {coverage.chain_residue_start}–{coverage.chain_residue_end} "
            f"of {coverage.full_protein_length} ({coverage.coverage_fraction:.0%} of full protein). "
            f"Regions: {', '.join(coverage.regions_covered)}."
        )

    if artifacts:
        interp_parts.append(f"Filtered {len(artifacts)} crystallization artifacts: {', '.join(artifacts[:5])}.")
    for s in sites:
        if s.num_residues > 0:
            interp_parts.append(f"Site {s.ligand_id}: {s.num_residues} contact residues (positions: {', '.join(str(p) for p in s.residue_positions[:10])}{'...' if len(s.residue_positions) > 10 else ''}).")
    site_types = set(s.site_type for s in sites)
    if "allosteric" in site_types:
        interp_parts.append("Allosteric site(s) detected — these may offer selectivity advantages.")
    if "cofactor" in site_types:
        interp_parts.append("Cofactor site(s) present — targeting these risks off-target effects on related enzymes.")

    result = BindingSiteMap(
        pdb_id=pdb_id.upper(),
        num_sites=len(sites),
        sites=sites,
        construct_coverage=coverage,
        artifact_ligands_filtered=artifacts,
        interpretation=" ".join(interp_parts),
    )
    return result.model_dump()


# ---------------------------------------------------------------------------
# Tool 4: get_ligand_history
# ---------------------------------------------------------------------------

@mcp.tool(
    name="GetLigandHistory",
    title="Get Ligand History",
    tags={"chemistry", "competitive-landscape"},
    annotations=READ_ONLY_ANNOTATIONS,
)
async def get_ligand_history(
    uniprot_id: str | None = None,
    pdb_id: str | None = None,
) -> dict:
    """Retrieve the bioactivity landscape for a drug target from ChEMBL.

    Shows what compounds have been tested, how potent the best ones are,
    whether any have reached clinical trials, and how crowded the
    competitive landscape is.

    Use this to decide whether to target KNOWN binding sites (where
    SAR exists) or seek NOVEL sites (where de novo design has an
    advantage). A crowded landscape suggests new modalities or
    allosteric approaches; an untargeted landscape suggests opportunity
    but less prior validation.

    Provide EITHER uniprot_id or pdb_id (UniProt preferred for accuracy).
    """
    if not uniprot_id and not pdb_id:
        return {"error": "Provide either uniprot_id or pdb_id"}

    _, uniprot_id = await _resolve_identifiers(pdb_id, uniprot_id)

    if not uniprot_id:
        return {"error": "Could not resolve UniProt ID"}

    # Map to ChEMBL
    target = await chembl.get_target_by_uniprot(uniprot_id)
    chembl_target_id = None
    if target:
        chembl_target_id = target.get("target_chembl_id")

    if not chembl_target_id:
        return LigandHistory(
            uniprot_id=uniprot_id,
            chembl_target_id=None,
            total_compounds_tested=0,
            total_assays=0,
            best_potency_nm=None,
            median_potency_nm=None,
            clinical_candidates=[],
            competitive_landscape="untargeted",
            interpretation=f"Target {uniprot_id} not found in ChEMBL. No bioactivity data available — this may be a novel target without pharmacological characterization.",
        ).model_dump()

    # Get bioactivities
    activities = await chembl.get_bioactivities(chembl_target_id)
    mechanisms = await chembl.get_clinical_candidates(chembl_target_id)
    analysis = analyze_bioactivities(activities, mechanisms)

    result = LigandHistory(
        uniprot_id=uniprot_id,
        chembl_target_id=chembl_target_id,
        **analysis,
    )
    return result.model_dump()


# ---------------------------------------------------------------------------
# Tool 5: check_conservation
# ---------------------------------------------------------------------------

# NCBI taxonomy IDs for preclinical model organisms
_ORGANISM_IDS: dict[str, int] = {
    "mouse": 10090,      # Mus musculus
    "rat": 10116,        # Rattus norvegicus
    "cynomolgus": 9541,  # Macaca fascicularis — closest primate preclinical model
}


@mcp.tool(
    name="CheckConservation",
    title="Check Conservation",
    tags={"conservation", "translatability"},
    annotations=READ_ONLY_ANNOTATIONS,
)
async def check_conservation(
    uniprot_id: str,
    residue_positions: list[int],
    species: list[str] | None = None,
) -> dict:
    """Check conservation at binding site residues across mouse, rat, and cynomolgus.

    Critical for preclinical translatability: if key binding site residues
    differ between human and a preclinical model, that species' efficacy
    data may not predict human response. Non-conserved positions are flagged
    for each species individually.

    Conservation > 90%: excellent — species should recapitulate human binding.
    Conservation 70-90%: acceptable — verify non-conserved positions are not
    critical contact residues.
    Conservation < 70%: caution for that species — consider a better-conserved
    alternative. Cynomolgus (macaque) is the closest primate model and often
    shows higher conservation than rodents when the target has primate-specific
    sequence features.

    Default species checked: mouse, rat, cynomolgus. Pass a custom `species`
    list to restrict or reorder (supported values: 'mouse', 'rat', 'cynomolgus').

    Provide the human UniProt accession and residue positions from GetBindingSites.
    """
    if not residue_positions:
        return {"error": "Provide a list of residue positions to check"}

    requested_species = species if species is not None else ["mouse", "rat", "cynomolgus"]

    # Fetch human sequence once
    try:
        human_seq = await uniprot.get_sequence(uniprot_id)
    except APIError as e:
        return {"error": f"Failed to fetch human sequence: {e}"}

    if not human_seq:
        return {"error": f"No sequence found for {uniprot_id}"}

    species_results: list[SpeciesConservation] = []

    for name in requested_species:
        org_id = _ORGANISM_IDS.get(name)
        if org_id is None:
            # Unknown species name — skip silently (or note it)
            species_results.append(SpeciesConservation(
                species=name,
                organism_id=0,
                note=f"Unknown species '{name}'. Supported: mouse, rat, cynomolgus.",
            ))
            continue

        # Resolve ortholog accession
        ortholog = await uniprot.get_ortholog(uniprot_id, org_id)
        if not ortholog:
            species_results.append(SpeciesConservation(
                species=name,
                organism_id=org_id,
                ortholog_uniprot=None,
                note=f"No ortholog found for {uniprot_id} in {name} (organism_id={org_id}).",
            ))
            continue

        # Fetch ortholog sequence
        try:
            ortholog_seq = await uniprot.get_sequence(ortholog)
        except APIError as e:
            species_results.append(SpeciesConservation(
                species=name,
                organism_id=org_id,
                ortholog_uniprot=ortholog,
                note=f"Found {name} ortholog {ortholog} but failed to retrieve sequence: {e}",
            ))
            continue

        # Per-position comparison using local context matching to handle indels
        conserved_count = 0
        non_conserved: list[ResidueConservation] = []

        for pos in residue_positions:
            idx = pos - 1  # Convert to 0-based
            human_res = human_seq[idx] if idx < len(human_seq) else "?"

            ortholog_res, map_score = _find_ortholog_residue(human_seq, ortholog_seq, idx)
            uncertain = map_score < 0.5

            is_same = human_res == ortholog_res
            is_conservative = _is_conservative_substitution(human_res, ortholog_res)

            if is_same:
                conserved_count += 1
            else:
                non_conserved.append(ResidueConservation(
                    position=pos,
                    human_residue=human_res,
                    ortholog_residue=ortholog_res,
                    is_conserved=False,
                    is_conservative_substitution=is_conservative,
                    mapping_uncertain=uncertain,
                ))

        fraction = conserved_count / len(residue_positions) if residue_positions else 0.0
        species_results.append(SpeciesConservation(
            species=name,
            organism_id=org_id,
            ortholog_uniprot=ortholog,
            residues_conserved=conserved_count,
            conservation_fraction=round(fraction, 3),
            non_conserved=non_conserved,
        ))

    # Build aggregate interpretation
    interp = _build_conservation_interpretation(uniprot_id, residue_positions, species_results)

    result = ConservationResult(
        human_uniprot=uniprot_id,
        residues_checked=len(residue_positions),
        species_results=species_results,
        interpretation=interp,
    )
    return result.model_dump()


def _build_conservation_interpretation(
    uniprot_id: str,
    residue_positions: list[int],
    species_results: list[SpeciesConservation],
) -> str:
    """Build a human-readable aggregate interpretation of multi-species conservation."""
    parts: list[str] = []

    # Per-species fractions summary
    frac_summaries = []
    for sr in species_results:
        if sr.ortholog_uniprot:
            pct = f"{sr.conservation_fraction:.0%}"
            frac_summaries.append(f"{sr.species} {pct}")
        else:
            frac_summaries.append(f"{sr.species} (no ortholog)")
    if frac_summaries:
        parts.append("Conservation at checked residues: " + ", ".join(frac_summaries) + ".")

    # Identify best model (highest fraction among species WITH an ortholog)
    with_ortholog = [sr for sr in species_results if sr.ortholog_uniprot]
    if with_ortholog:
        best = max(with_ortholog, key=lambda sr: sr.conservation_fraction)
        parts.append(
            f"Recommended preclinical model: {best.species} "
            f"({best.conservation_fraction:.0%} conservation, UniProt {best.ortholog_uniprot})."
        )

    # Flag any species below 0.7
    low = [sr for sr in with_ortholog if sr.conservation_fraction < 0.7]
    for sr in low:
        non_pos = [str(r.position) for r in sr.non_conserved]
        parts.append(
            f"Caution: {sr.species} conservation is low ({sr.conservation_fraction:.0%}). "
            f"Non-conserved positions: {', '.join(non_pos) if non_pos else 'none listed'}. "
            f"{sr.species.capitalize()} models may not recapitulate human target binding at these residues."
        )

    # Flag species without orthologs
    missing = [sr for sr in species_results if not sr.ortholog_uniprot]
    for sr in missing:
        parts.append(f"No ortholog found for {uniprot_id} in {sr.species} — conservation cannot be assessed.")

    # Flag uncertain mappings across all species
    total_uncertain = sum(
        sum(1 for r in sr.non_conserved if r.mapping_uncertain)
        for sr in with_ortholog
    )
    if total_uncertain > 0:
        parts.append(
            f"{total_uncertain} position(s) across species had low-confidence ortholog mapping — treat those calls cautiously."
        )

    return " ".join(parts) if parts else "No conservation data available."


# ---------------------------------------------------------------------------
# Tool 6: check_known_variants
# ---------------------------------------------------------------------------

@mcp.tool(
    name="CheckKnownVariants",
    title="Check Known Variants",
    tags={"variants", "resistance"},
    annotations=READ_ONLY_ANNOTATIONS,
)
async def check_known_variants(uniprot_id: str, residue_positions: list[int]) -> dict:
    """Flag known sequence variants and mutagenesis hits at binding-site residues.

    Binding-site residues that are documented disease/resistance variants (e.g.
    EGFR T790M) mark pockets that mutate under drug pressure — a key risk signal
    when choosing where to design. Source: UniProt Natural variant + Mutagenesis
    features. Provide the human UniProt accession and positions from GetBindingSites.
    """
    if not residue_positions:
        return {"error": "Provide a list of residue positions to check"}
    try:
        entry = await uniprot.get_entry(uniprot_id)
    except APIError as e:
        return {"error": f"Failed to fetch UniProt entry: {e}"}
    found = parse_known_variants(entry, set(residue_positions))
    variants = [KnownVariant(**v) for v in found]
    if variants:
        interp = (f"{len(variants)} known variant(s) at the queried binding-site residues: "
                  + "; ".join(f"{v.original_residue}{v.position}{v.variant_residue} ({v.feature_type})" for v in variants)
                  + ". Variants at contact residues can drive resistance or alter binding — design with this in mind.")
    else:
        interp = "No known sequence variants or mutagenesis records at the queried residues — the pocket appears genetically stable in UniProt."
    return VariantCheckResult(
        uniprot_id=uniprot_id,
        positions_checked=len(residue_positions),
        variants=variants,
        interpretation=interp,
    ).model_dump()


# ---------------------------------------------------------------------------
# Tool 7: search_target_literature
# ---------------------------------------------------------------------------

@mcp.tool(
    name="SearchTargetLiterature",
    title="Search Target Literature",
    tags={"literature"},
    annotations=READ_ONLY_ANNOTATIONS,
)
async def search_target_literature(
    gene_name: str,
    context: str | None = None,
    max_results: int = 10,
) -> dict:
    """Search PubMed for recent structural biology and drug design papers.

    Focuses specifically on binding site characterization, allosteric
    mechanisms, resistance mutations, and prior design campaigns — the
    literature most relevant to planning a new binder design effort.

    Use the optional `context` parameter to narrow results, e.g.:
    - context='allosteric' for allosteric site literature
    - context='resistance' for resistance mutation papers
    - context='antibody' for biologic-focused papers
    - context='oncology' for disease-specific context

    Call this LAST — after CharacterizeTarget, GetRelatedStructures,
    GetBindingSites, GetLigandHistory, and CheckConservation — to see if
    the literature reveals insights not captured in database records (e.g.,
    cryptic sites found by MD simulation, unpublished allosteric mechanisms).
    """
    try:
        result = await pubmed.search_target_literature(
            gene_name=gene_name,
            context=context,
            max_results=max_results,
        )
    except APIError as e:
        return {"error": f"PubMed search failed: {e}"}

    papers = [
        PaperResult(**p)
        for p in result.get("papers", [])
    ]

    total = result.get("total_found", 0)
    interp = f"Found {total} papers matching structural/design criteria for {gene_name}."
    if total > 100:
        interp += " This is a well-studied target with extensive structural literature."
    elif total < 10:
        interp += " Limited structural literature — this target may be under-explored or recently discovered."

    lit_result = LiteratureResult(
        query=result.get("query", ""),
        total_found=total,
        papers=papers,
        interpretation=interp,
    )
    return lit_result.model_dump()


# ---------------------------------------------------------------------------
# Clustering helper (pure function — unit-tested independently)
# ---------------------------------------------------------------------------

def cluster_pockets(sites: list[dict], jaccard_threshold: float = 0.3) -> list[dict]:
    """Greedy single-link clustering of pockets across structures by residue overlap.

    Each site dict must have: pdb_id, ligand_id, site_type, residue_positions.
    Sites with no residue_positions are skipped.
    Returns clusters sorted by structure_count descending.
    """
    clusters: list[dict] = []
    for site in sites:
        rset = set(site["residue_positions"])
        if not rset:
            continue
        placed = False
        for cl in clusters:
            inter = len(cl["_residues"] & rset)
            union = len(cl["_residues"] | rset)
            if union and inter / union >= jaccard_threshold:
                cl["_residues"] |= rset
                cl["occurrences"].append((site["pdb_id"], site["ligand_id"]))
                cl["_types"].append(site["site_type"])
                placed = True
                break
        if not placed:
            clusters.append({
                "_residues": set(rset),
                "occurrences": [(site["pdb_id"], site["ligand_id"])],
                "_types": [site["site_type"]],
            })
    out = []
    for cl in clusters:
        types = cl["_types"]
        out.append({
            "residue_union": sorted(cl["_residues"]),
            "occurrences": cl["occurrences"],
            "structure_count": len({pid for pid, _ in cl["occurrences"]}),
            "representative_site_type": max(set(types), key=types.count),
        })
    out.sort(key=lambda c: c["structure_count"], reverse=True)
    return out


# ---------------------------------------------------------------------------
# Tool 8: consolidate_binding_sites
# ---------------------------------------------------------------------------

@mcp.tool(name="ConsolidateBindingSites", title="Consolidate Binding Sites",
          tags={"structure", "binding-site"}, annotations=READ_ONLY_ANNOTATIONS)
async def consolidate_binding_sites(
    uniprot_id: str | None = None,
    pdb_id: str | None = None,
    limit: int = 10,
) -> dict:
    """Map the union of binding pockets across all structures of a target.

    Fans out GetBindingSites over the top structures and clusters pockets by
    residue overlap, so recurrent pockets (e.g. the ATP site appearing in most
    structures) stand out from one-off or artifact sites. The heaviest tool —
    downloads several coordinate files. Provide uniprot_id (preferred) or pdb_id.
    """
    if not uniprot_id and not pdb_id:
        return {"error": "Provide either uniprot_id or pdb_id"}

    _, uniprot_id = await _resolve_identifiers(pdb_id, uniprot_id)

    if not uniprot_id:
        return {"error": "Could not resolve UniProt ID"}

    try:
        pdb_ids = await pdb.search_by_uniprot(uniprot_id, limit=limit)
    except APIError as e:
        return {"error": f"PDB search failed: {e}"}

    import asyncio
    sem = asyncio.Semaphore(6)

    async def _sites(pid: str) -> list[dict]:
        async with sem:
            try:
                data = await pdb.get_binding_sites(pid)
            except APIError:
                return []
        out = []
        for lig in data.get("ligands", []):
            positions = sorted({
                c["residue_number"]
                for c in lig.get("contacts", [])
                if c["residue_number"] > 0
            })
            if positions:
                out.append({
                    "pdb_id": pid.upper(),
                    "ligand_id": lig.get("comp_id", ""),
                    "site_type": _classify_site_type(
                        lig.get("comp_id", ""), lig.get("name", "")
                    ),
                    "residue_positions": positions,
                })
        return out

    per_structure = await asyncio.gather(*[_sites(p) for p in pdb_ids])
    all_sites = [s for group in per_structure for s in group]
    clusters = cluster_pockets(all_sites)
    pockets = [ConsolidatedPocket(**c) for c in clusters]

    if pockets:
        top = pockets[0]
        interp = (
            f"{len(pockets)} distinct pocket(s) across {len(pdb_ids)} structures. "
            f"Most recurrent: a {top.representative_site_type} pocket in {top.structure_count} structure(s) "
            f"(residues {top.residue_union[0]}–{top.residue_union[-1]}) — the dominant, best-validated site."
        )
    else:
        interp = f"No non-artifact pockets found across {len(pdb_ids)} structures."

    return ConsolidationResult(
        uniprot_id=uniprot_id,
        structures_analyzed=len(pdb_ids),
        pockets=pockets,
        interpretation=interp,
    ).model_dump()


# ---------------------------------------------------------------------------
# Orchestration prompt
# ---------------------------------------------------------------------------

@mcp.prompt()
def binding_site_assessment(
    pdb_id: str,
    indication: str | None = None,
) -> str:
    """Guide a systematic binding site evaluation for de novo binder design.

    This prompt orchestrates all PocketScout tools in scientific workflow
    order to produce a ranked assessment of candidate binding regions.

    The output should include:
    - Target characterization with druggability context
    - Complete binding site map across available structures
    - Competitive landscape assessment
    - Cross-species conservation for translatability
    - Literature-informed insights
    - Ranked recommendation of binding regions for de novo design
    """
    context_line = f" in the context of {indication}" if indication else ""

    return f"""Evaluate PDB structure {pdb_id} as a target for de novo binder design{context_line}.

Follow this systematic workflow using PocketScout tools:

**Step 1 — Target Context**
Call CharacterizeTarget with pdb_id="{pdb_id}".
Establish: What is this protein? What family? Where is it located?
Note any AlphaFold confidence warnings for later cross-referencing.

**Step 2 — Structural Coverage**
Call GetRelatedStructures with pdb_id="{pdb_id}".
Assess: How many structures exist? What ligands? What quality?
Identify the best structures for binding site analysis.

**Step 3 — Binding Site Map**
Call GetBindingSites on the most informative structure(s).
For well-studied targets, analyze 2-3 structures with different ligands
to capture multiple pocket conformations.
Map: Known binding sites, classification, druggability.

**Step 4 — Competitive Landscape**
Call GetLigandHistory using the UniProt ID from Step 1.
Assess: Is this target crowded or greenfield?
What modalities have been tried? Where is the opportunity?

**Step 5 — Translatability**
Call CheckConservation with binding site residue positions from Step 3.
Flag: Are key contact residues conserved in mouse?
What are the implications for preclinical models?

**Step 6 — Literature Context**
Call SearchTargetLiterature with the gene name from Step 1.
Look for: Cryptic/allosteric sites, resistance mutations, recent
structural insights not captured in database records.

**Synthesis — Ranked Recommendation**
Combine all evidence into a ranked assessment:
1. Which binding region(s) are most promising for de novo design?
2. What are the trade-offs between sites (druggability, selectivity,
   conservation, competitive landscape)?
3. What modality is recommended for each site (small molecule,
   biologic, de novo protein binder)?
4. What parameters would you recommend for downstream computational
   design (target residues, constraints, flexibility considerations)?
5. What are the key risks, unknowns, or additional experiments needed?

Be specific and evidence-based. Cite PDB IDs, ChEMBL data, and
literature where relevant. Flag uncertainties clearly."""


# ---------------------------------------------------------------------------
# Quick triage prompt
# ---------------------------------------------------------------------------

@mcp.prompt()
def target_briefing(target: str) -> str:
    """Quick triage briefing for a drug target — what it is and where its pockets are."""
    return f"""Give me a fast binding-site briefing on {target} using PocketScout.

1. CharacterizeTarget — what is this protein (family, location, structure confidence)?
2. GetRelatedStructures — how well-characterized is it structurally?
3. GetBindingSites on the best structure — what pockets exist and how druggable are they?
4. GetLigandHistory — is the chemical space crowded or open?
5. SearchTargetLiterature — any notable allosteric/cryptic sites or resistance themes?

Summarize in a short briefing: what the target is, the main pocket(s), the
competitive landscape, and the one or two things worth knowing before going
deeper. Keep it concise — this is triage, not a full design workup."""


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

async def _compute_construct_coverage(pdb_id: str) -> ConstructCoverage | None:
    """Determine what portion of the full protein a PDB structure covers.

    Cross-references the PDB chain residue range against UniProt topology
    to determine whether the structure covers extracellular, transmembrane,
    or cytoplasmic regions. This is the key check for de novo binder
    design feasibility.
    """
    pdb_id = pdb_id.upper()

    # Get the UniProt residue range from PDB
    residue_range = await pdb.get_uniprot_residue_range(pdb_id)
    if not residue_range:
        return None

    chain_start, chain_end = residue_range

    # Get UniProt mapping to fetch topology
    mappings = await pdb.get_uniprot_mapping(pdb_id)
    if not mappings:
        return None

    uniprot_id = mappings[0]

    try:
        entry = await uniprot.get_entry(uniprot_id)
    except APIError:
        return None

    profile_data = parse_target_profile(entry)
    seq_length = profile_data["sequence_length"]
    topology = profile_data.get("topology", [])

    coverage_frac = (chain_end - chain_start + 1) / seq_length if seq_length > 0 else 0.0

    # Determine which topological regions the construct covers
    regions_covered = set()
    if not topology:
        # No topology annotations — protein is likely entirely soluble/cytoplasmic
        regions_covered.add("cytoplasmic (no membrane topology annotated)")
    else:
        for topo in topology:
            topo_start = topo["start"]
            topo_end = topo["end"]
            region_type = topo["region_type"]

            if region_type == "signal_peptide":
                continue

            # Check for overlap between construct and topology region
            overlap_start = max(chain_start, topo_start)
            overlap_end = min(chain_end, topo_end)
            if overlap_start <= overlap_end:
                regions_covered.add(region_type)

    # Generate accessibility warning
    warning = ""
    regions_list = sorted(regions_covered)

    if topology and "extracellular" not in regions_covered:
        if "cytoplasmic" in regions_covered or any("cytoplasmic" in r for r in regions_covered):
            warning = (
                f"This structure covers ONLY intracellular/cytoplasmic regions (residues {chain_start}–{chain_end}). "
                f"De novo protein binders CANNOT access intracellular targets — they act extracellularly. "
                f"For protein binder design, use structures covering the extracellular domain instead. "
                f"This structure is suitable for small molecule or degrader design only."
            )
        elif "transmembrane" in regions_covered:
            warning = (
                f"This structure covers only the transmembrane region. "
                f"De novo protein binder design requires extracellular domain structures."
            )

    return ConstructCoverage(
        pdb_id=pdb_id,
        chain_residue_start=chain_start,
        chain_residue_end=chain_end,
        full_protein_length=seq_length,
        coverage_fraction=round(coverage_frac, 3),
        regions_covered=regions_list,
        accessibility_warning=warning,
    )


def _find_ortholog_residue(
    human_seq: str, ortholog_seq: str, human_idx: int, window: int = 7
) -> tuple[str, float]:
    """Find the ortholog residue corresponding to a human sequence position.

    Uses local context matching: extracts a window around the human position,
    searches for the best match in the ortholog sequence, and returns the
    corresponding central residue. This handles insertions/deletions that
    shift numbering between orthologs.

    Returns a tuple of (residue, match_fraction) where match_fraction is
    the fraction of the context window that matched at the best position.
    A low match_fraction (<0.5) indicates an uncertain mapping.
    """
    if human_idx >= len(human_seq):
        return "?", 0.0

    # Extract context window from human sequence
    start = max(0, human_idx - window)
    end = min(len(human_seq), human_idx + window + 1)
    human_context = human_seq[start:end]
    offset_in_context = human_idx - start

    # Guard against empty context (degenerate case)
    if len(human_context) == 0:
        return "?", 0.0

    # Search for best matching position in ortholog sequence
    best_score = -1
    best_ortholog_idx = human_idx  # fallback to same position

    # Search within a reasonable range around the expected position
    search_start = max(0, human_idx - 50)
    search_end = min(len(ortholog_seq) - len(human_context) + 1, human_idx + 50)

    for i in range(search_start, search_end):
        ortholog_window = ortholog_seq[i:i + len(human_context)]
        if len(ortholog_window) < len(human_context):
            continue
        score = sum(1 for a, b in zip(human_context, ortholog_window) if a == b)
        if score > best_score:
            best_score = score
            best_ortholog_idx = i + offset_in_context

    match_fraction = round(best_score / len(human_context), 3) if best_score >= 0 else 0.0

    if best_ortholog_idx < len(ortholog_seq):
        return ortholog_seq[best_ortholog_idx], match_fraction
    return "?", match_fraction


# Conservative amino acid substitution groups
_CONSERVATIVE_GROUPS = [
    frozenset("MILV"),       # Hydrophobic aliphatic
    frozenset("FYW"),        # Aromatic
    frozenset("KRH"),        # Positive charged
    frozenset("DE"),         # Negative charged
    frozenset("STNQ"),       # Polar uncharged
    frozenset("AG"),         # Small
]


def _is_conservative_substitution(aa1: str, aa2: str) -> bool:
    """Check if two amino acids are conservative substitutions."""
    if aa1 == aa2:
        return True
    for group in _CONSERVATIVE_GROUPS:
        if aa1 in group and aa2 in group:
            return True
    return False


# Known cofactors and allosteric modulators
_COFACTOR_IDS = {"ATP", "ADP", "AMP", "NAD", "NAP", "FAD", "FMN", "SAM", "COA", "HEM", "GTP", "GDP"}
_KNOWN_ALLOSTERIC_KEYWORDS = {"allosteric", "regulatory", "exosite"}


def _classify_site_type(comp_id: str, ligand_name: str) -> str:
    """Classify a binding site based on the ligand identity."""
    if comp_id in _COFACTOR_IDS:
        return "cofactor"
    name_lower = ligand_name.lower()
    for kw in _KNOWN_ALLOSTERIC_KEYWORDS:
        if kw in name_lower:
            return "allosteric"
    # Default assumption: drug-like ligands are at the orthosteric site
    return "orthosteric"


def _assess_druggability(site_type: str, comp_id: str, num_residues: int = 0) -> str:
    """Generate druggability notes based on site classification and pocket size."""
    size_note = ""
    if num_residues > 0:
        if num_residues < 12:
            size_note = f" Compact pocket ({num_residues} contact residues) — suitable for small molecule design."
        elif num_residues <= 25:
            size_note = f" Medium pocket ({num_residues} contact residues) — suitable for small molecules or peptide mimetics."
        else:
            size_note = f" Large binding interface ({num_residues} contact residues) — consider de novo protein binders or biologics."

    if site_type == "cofactor":
        return (
            f"Cofactor binding site ({comp_id}). Targeting cofactor sites risks selectivity issues "
            "as many related enzymes share similar cofactor-binding pockets. Consider allosteric "
            f"or substrate-competitive approaches instead.{size_note}"
        )
    elif site_type == "allosteric":
        return (
            "Allosteric site — potentially high selectivity since allosteric pockets are less "
            "conserved across protein families than orthosteric sites. Good candidate for "
            f"de novo design. May require experimental validation of functional effect.{size_note}"
        )
    else:
        return (
            f"Orthosteric/primary binding site (ligand: {comp_id}). Well-characterized pocket.{size_note}"
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """Run the PocketScout MCP server.

    Supports two transport modes:
      - stdio (default): for local MCP clients like Claude Desktop
      - http: serves over streamable HTTP for remote access

    Usage:
      pocketscout-mcp          # stdio mode
      pocketscout-mcp http     # HTTP mode on port 8000
    """
    import sys
    import os

    transport = sys.argv[1] if len(sys.argv) > 1 else "stdio"

    if transport == "http":
        host = os.environ.get("HOST", "0.0.0.0")
        port = int(os.environ.get("PORT", "8000"))
        mcp.run(transport="streamable-http", host=host, port=port)
    else:
        mcp.run(transport="stdio")


if __name__ == "__main__":
    main()

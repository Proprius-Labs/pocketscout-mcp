"""PocketScout MCP Server

Scout the binding landscape before you design the binder.

An MCP server that aggregates structural, chemical, and literature data
to evaluate druggable pockets on protein targets before computational
design. Designed to fill the gap between "I have a target" and
"I'm running RFdiffusion."

Tools are ordered to reflect the natural scientific workflow:
1. characterize_target — biological context + AlphaFold confidence
2. get_related_structures — structural coverage for the target
3. get_binding_sites — map known pockets from co-crystals
4. get_ligand_history — competitive landscape from ChEMBL
5. check_conservation — human vs. mouse translatability
6. search_target_literature — recent structural/design insights

The binding_site_assessment prompt orchestrates all tools into a
ranked recommendation.
"""

from __future__ import annotations

from fastmcp import FastMCP

from .clients.uniprot import UniProtClient, parse_target_profile
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
    ConservationResult,
    PaperResult,
    LiteratureResult,
)

# ---------------------------------------------------------------------------
# Server + clients
# ---------------------------------------------------------------------------

mcp = FastMCP(
    "PocketScout",
    instructions=(
        "Binding site intelligence for drug target assessment. "
        "Aggregates structural, chemical, and literature data to evaluate "
        "druggable pockets before computational design. "
        "Use the binding_site_assessment prompt for a complete workflow."
    ),
)

uniprot = UniProtClient()
pdb = PDBClient()
alphafold = AlphaFoldClient()
chembl = ChEMBLClient()
pubmed = PubMedClient()


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

@mcp.tool
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
        af_data = analyze_confidence(af_prediction, profile_data["sequence_length"])
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

@mcp.tool
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

    Call this AFTER characterize_target and BEFORE get_binding_sites to
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

    if not pdb_ids:
        return {
            "uniprot_id": uniprot_id,
            "gene_name": "",
            "total_structures": 0,
            "structures": [],
            "unique_ligands": [],
            "has_apo_structure": False,
            "best_resolution": None,
            "interpretation": "No experimental structures found in PDB. Consider using AlphaFold predicted structure, but note that predicted binding sites have lower reliability than those from experimental structures.",
        }

    # Fetch metadata for each structure (limit to top entries to avoid slow response)
    structures = []
    all_ligands = set()
    has_apo = False

    for pid in pdb_ids[:limit]:
        try:
            entry = await pdb.get_entry(pid)
            meta = parse_structure_metadata(entry)

            # Get ligands
            np_entities = await pdb.get_nonpolymer_entities(pid)
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
                    all_ligands.add(comp_id)

            if not ligand_ids:
                has_apo = True

            organism = "unknown"
            entities = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])

            structures.append(RelatedStructure(
                pdb_id=pid,
                title=meta["title"],
                resolution=meta["resolution"],
                method=meta["method"],
                organism=organism,
                ligand_ids=ligand_ids,
                release_date=meta["release_date"],
            ))
        except APIError:
            continue

    # Sort by resolution
    structures.sort(key=lambda s: s.resolution if s.resolution else 99.0)
    best_res = structures[0].resolution if structures else None

    # Get gene name from UniProt (quick)
    gene_name = ""
    try:
        entry = await uniprot.get_entry(uniprot_id)
        gene_name = parse_target_profile(entry)["gene_name"]
    except APIError:
        pass

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

@mcp.tool
async def get_binding_sites(pdb_id: str) -> dict:
    """Map all known binding sites in a protein structure from co-crystallized ligands.

    Identifies binding pockets by analyzing non-polymer entities (ligands,
    cofactors) in the structure, filtering out crystallization artifacts
    (glycerol, PEG, sulfate, etc.), and classifying each site by type.

    Each site includes druggability assessment and modality recommendations:
    - Small deep pockets (< 20 contact residues) favor small molecules
    - Large flat interfaces (> 30 residues) favor biologics or de novo protein binders
    - Allosteric sites may offer selectivity advantages over orthosteric sites

    Call this on specific PDB IDs identified by get_related_structures.
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
            interpretation=f"No non-artifact ligands found in {pdb_id}. This is an apo structure or contains only crystallization artifacts ({', '.join(artifacts) if artifacts else 'none detected'}). Check other structures for this target via get_related_structures.",
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

@mcp.tool
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

@mcp.tool
async def check_conservation(
    uniprot_id: str,
    residue_positions: list[int],
) -> dict:
    """Check human vs. mouse conservation at specific binding site residues.

    Critical for preclinical translatability: if key binding site residues
    differ between human and mouse, mouse efficacy models may not predict
    human response. Non-conserved positions are specifically flagged.

    Conservation > 90%: excellent — mouse model should recapitulate binding.
    Conservation 70-90%: acceptable — verify non-conserved positions aren't
    in the binding interface.
    Conservation < 70%: caution — consider rat, cyno, or humanized models.

    Provide the human UniProt accession and a list of residue positions
    (from get_binding_sites) to check.
    """
    if not residue_positions:
        return {"error": "Provide a list of residue positions to check"}

    # Get human sequence
    try:
        human_seq = await uniprot.get_sequence(uniprot_id)
    except APIError as e:
        return {"error": f"Failed to fetch human sequence: {e}"}

    if not human_seq:
        return {"error": f"No sequence found for {uniprot_id}"}

    # Find mouse ortholog
    mouse_id = await uniprot.get_mouse_ortholog(uniprot_id)
    if not mouse_id:
        return ConservationResult(
            human_uniprot=uniprot_id,
            mouse_uniprot=None,
            residues_checked=len(residue_positions),
            residues_conserved=0,
            conservation_fraction=0.0,
            non_conserved=[],
            interpretation=f"No mouse ortholog found for {uniprot_id}. Cannot assess cross-species conservation. This may indicate the target is not conserved in rodents, or the ortholog has a different gene name.",
        ).model_dump()

    try:
        mouse_seq = await uniprot.get_sequence(mouse_id)
    except APIError:
        return ConservationResult(
            human_uniprot=uniprot_id,
            mouse_uniprot=mouse_id,
            residues_checked=len(residue_positions),
            residues_conserved=0,
            conservation_fraction=0.0,
            non_conserved=[],
            interpretation=f"Found mouse ortholog {mouse_id} but failed to retrieve sequence.",
        ).model_dump()

    # Compare residues using local context matching to handle insertions/deletions.
    # For each human position, we extract a window of surrounding sequence and
    # find where it occurs in the mouse sequence, then compare the central residue.
    conserved_count = 0
    non_conserved = []

    for pos in residue_positions:
        idx = pos - 1  # Convert to 0-based
        human_res = human_seq[idx] if idx < len(human_seq) else "?"

        # Find the corresponding mouse residue by context matching
        mouse_res = _find_ortholog_residue(human_seq, mouse_seq, idx)

        is_same = human_res == mouse_res
        is_conservative = _is_conservative_substitution(human_res, mouse_res)

        if is_same:
            conserved_count += 1
        else:
            non_conserved.append(ResidueConservation(
                position=pos,
                human_residue=human_res,
                mouse_residue=mouse_res,
                is_conserved=False,
                is_conservative_substitution=is_conservative,
            ))

    fraction = conserved_count / len(residue_positions) if residue_positions else 0.0

    # Generate interpretation
    if fraction > 0.9:
        interp = f"Excellent conservation ({fraction:.0%}). Mouse models should faithfully recapitulate human target binding."
    elif fraction > 0.7:
        non_cons_positions = [str(r.position) for r in non_conserved]
        interp = (
            f"Good conservation ({fraction:.0%}). "
            f"Non-conserved positions at: {', '.join(non_cons_positions)}. "
            "Verify these positions are not critical contact residues in your binding site of interest."
        )
    else:
        interp = (
            f"Low conservation ({fraction:.0%}). "
            "Significant differences between human and mouse at binding site residues. "
            "Mouse efficacy models may not predict human response. Consider cynomolgus monkey, "
            "humanized mouse models, or other species with better conservation."
        )

    result = ConservationResult(
        human_uniprot=uniprot_id,
        mouse_uniprot=mouse_id,
        residues_checked=len(residue_positions),
        residues_conserved=conserved_count,
        conservation_fraction=round(fraction, 3),
        non_conserved=non_conserved,
        interpretation=interp,
    )
    return result.model_dump()


# ---------------------------------------------------------------------------
# Tool 6: search_target_literature
# ---------------------------------------------------------------------------

@mcp.tool
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

    Call this LAST — after structural and chemical data — to see if the
    literature reveals insights not captured in database records (e.g.,
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
Call characterize_target with pdb_id="{pdb_id}".
Establish: What is this protein? What family? Where is it located?
Note any AlphaFold confidence warnings for later cross-referencing.

**Step 2 — Structural Coverage**
Call get_related_structures with pdb_id="{pdb_id}".
Assess: How many structures exist? What ligands? What quality?
Identify the best structures for binding site analysis.

**Step 3 — Binding Site Map**
Call get_binding_sites on the most informative structure(s).
For well-studied targets, analyze 2-3 structures with different ligands
to capture multiple pocket conformations.
Map: Known binding sites, classification, druggability.

**Step 4 — Competitive Landscape**
Call get_ligand_history using the UniProt ID from Step 1.
Assess: Is this target crowded or greenfield?
What modalities have been tried? Where is the opportunity?

**Step 5 — Translatability**
Call check_conservation with binding site residue positions from Step 3.
Flag: Are key contact residues conserved in mouse?
What are the implications for preclinical models?

**Step 6 — Literature Context**
Call search_target_literature with the gene name from Step 1.
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
    human_seq: str, mouse_seq: str, human_idx: int, window: int = 7
) -> str:
    """Find the mouse residue corresponding to a human sequence position.

    Uses local context matching: extracts a window around the human position,
    searches for the best match in the mouse sequence, and returns the
    corresponding central residue. This handles insertions/deletions that
    shift numbering between orthologs.
    """
    if human_idx >= len(human_seq):
        return "?"

    # Extract context window from human sequence
    start = max(0, human_idx - window)
    end = min(len(human_seq), human_idx + window + 1)
    human_context = human_seq[start:end]
    offset_in_context = human_idx - start

    # Search for best matching position in mouse sequence
    best_score = -1
    best_mouse_idx = human_idx  # fallback to same position

    # Search within a reasonable range around the expected position
    search_start = max(0, human_idx - 50)
    search_end = min(len(mouse_seq) - len(human_context) + 1, human_idx + 50)

    for i in range(search_start, search_end):
        mouse_window = mouse_seq[i:i + len(human_context)]
        if len(mouse_window) < len(human_context):
            continue
        score = sum(1 for a, b in zip(human_context, mouse_window) if a == b)
        if score > best_score:
            best_score = score
            best_mouse_idx = i + offset_in_context

    if best_mouse_idx < len(mouse_seq):
        return mouse_seq[best_mouse_idx]
    return "?"


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
    """Run the PocketScout MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()

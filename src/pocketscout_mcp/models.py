"""Structured data models for PocketScout tool returns.

Field descriptions propagate into MCP schemas — Claude reads these
when deciding how to interpret tool results. Every description should
encode scientific reasoning, not just data labels.
"""

from __future__ import annotations

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Tool 1: characterize_target
# ---------------------------------------------------------------------------

class ConfidenceRegion(BaseModel):
    """A contiguous region of the protein with a confidence assessment."""
    start: int
    end: int
    mean_plddt: float = Field(description="Mean AlphaFold pLDDT score for this region. >90 = high confidence, 70-90 = moderate, <70 = low — structure predictions in low-confidence regions should not be trusted for binding site analysis.")
    assessment: str = Field(description="'high' (>90), 'moderate' (70-90), or 'low' (<70)")


class TargetProfile(BaseModel):
    """Biological context for a drug target protein."""
    uniprot_id: str
    gene_name: str
    protein_name: str
    organism: str
    function_summary: str = Field(description="Brief functional description from UniProt")
    protein_family: str = Field(description="Protein family classification, e.g. 'Protein kinase', 'GPCR family 1', 'Immunoglobulin superfamily'")
    subcellular_location: str = Field(description="Where the protein resides — critical for modality selection. Extracellular/membrane proteins are accessible to biologics; intracellular targets require small molecules or cell-penetrating designs.")
    disease_associations: list[str] = Field(default_factory=list, description="Diseases linked to this target in UniProt")
    target_class: str = Field(description="Druggability class: 'kinase', 'gpcr', 'protease', 'ion_channel', 'nuclear_receptor', 'ppi', 'enzyme_other', or 'other'")
    sequence_length: int = Field(description="Full protein sequence length in amino acids")
    alphafold_confidence: list[ConfidenceRegion] = Field(
        default_factory=list,
        description="AlphaFold pLDDT confidence by region. Low-confidence regions (<70) likely represent disordered or flexible segments where predicted binding sites should not be trusted."
    )
    alphafold_overall_confidence: float | None = Field(default=None, description="Global mean pLDDT across the full structure")
    low_confidence_warnings: list[str] = Field(
        default_factory=list,
        description="Specific warnings about regions where AlphaFold predictions are unreliable"
    )


# ---------------------------------------------------------------------------
# Tool 2: get_binding_sites
# ---------------------------------------------------------------------------

class BindingSiteResidue(BaseModel):
    """A residue involved in a binding interaction."""
    chain: str
    residue_name: str = Field(description="Three-letter amino acid code, e.g. 'ALA', 'GLU'")
    residue_number: int
    contact_type: str = Field(default="unknown", description="Type of interaction: 'hydrogen_bond', 'hydrophobic', 'ionic', 'pi_stacking', or 'unknown'")


class BindingSite(BaseModel):
    """A characterized binding site from a co-crystal structure."""
    site_id: str
    residues: list[BindingSiteResidue] = Field(description="Residues lining or contacting this pocket")
    residue_positions: list[int] = Field(description="Residue numbers only — for cross-referencing with conservation and confidence data")
    ligand_id: str | None = Field(default=None, description="PDB ligand ID (3-letter code) bound at this site, e.g. 'STI' for imatinib")
    ligand_name: str | None = Field(default=None, description="Common name of the ligand if known")
    site_type: str = Field(description="'orthosteric' (primary/active site), 'allosteric' (regulatory), 'cofactor' (catalytic cofactor like ATP/NAD), 'interface' (protein-protein), or 'unknown'")
    evidence: str = Field(default="co-crystal", description="How this site was identified: 'co-crystal', 'competition_assay', 'mutagenesis', 'computational'")
    num_residues: int = Field(description="Number of residues lining the pocket — rough proxy for pocket size")
    druggability_notes: str = Field(description="Assessment of suitability for different modalities. Small deep pockets favor small molecules; large flat interfaces favor biologics or de novo protein design.")


class BindingSiteMap(BaseModel):
    """Complete binding site map for a structure."""
    pdb_id: str
    num_sites: int
    sites: list[BindingSite]
    artifact_ligands_filtered: list[str] = Field(
        default_factory=list,
        description="Ligands excluded as crystallization artifacts (glycerol, PEG, sulfate, etc.)"
    )
    interpretation: str = Field(description="Summary of the binding landscape — how many distinct sites, which are most druggable, notable features")


# ---------------------------------------------------------------------------
# Tool 3: get_related_structures
# ---------------------------------------------------------------------------

class RelatedStructure(BaseModel):
    """A PDB structure related to the target protein."""
    pdb_id: str
    title: str
    resolution: float | None = Field(default=None, description="Resolution in Angstroms. <2.0 = high quality, 2.0-3.0 = acceptable, >3.0 = use with caution for binding site analysis.")
    method: str = Field(description="Experimental method: 'X-RAY DIFFRACTION', 'ELECTRON MICROSCOPY', 'SOLUTION NMR'")
    organism: str = Field(default="unknown")
    ligand_ids: list[str] = Field(default_factory=list, description="Ligands co-crystallized in this structure")
    release_date: str | None = Field(default=None)


class RelatedStructuresResult(BaseModel):
    """All available structures for a target protein."""
    uniprot_id: str
    gene_name: str
    total_structures: int
    structures: list[RelatedStructure] = Field(description="Structures sorted by resolution (best first)")
    unique_ligands: list[str] = Field(default_factory=list, description="Distinct ligands across all structures")
    has_apo_structure: bool = Field(description="Whether an unliganded structure exists — useful for comparing bound vs. unbound conformations")
    best_resolution: float | None = Field(default=None, description="Best available resolution across all structures")
    interpretation: str = Field(description="Summary: structural coverage, ligand diversity, quality assessment")


# ---------------------------------------------------------------------------
# Tool 4: get_ligand_history
# ---------------------------------------------------------------------------

class LigandHistory(BaseModel):
    """Bioactivity landscape from ChEMBL for a target."""
    uniprot_id: str
    chembl_target_id: str | None = Field(default=None)
    total_compounds_tested: int
    total_assays: int
    best_potency_nm: float | None = Field(default=None, description="Best reported potency (IC50/Ki/Kd) in nanomolar. <10 nM = very potent, 10-100 nM = potent, 100-1000 nM = moderate.")
    median_potency_nm: float | None = Field(default=None, description="Median potency across all tested compounds")
    clinical_candidates: list[str] = Field(default_factory=list, description="Compounds that have reached clinical trials")
    competitive_landscape: str = Field(
        description="Assessment: 'crowded' (>100 compounds with <1μM activity), "
                    "'moderate' (10-100), 'emerging' (<10 active compounds), "
                    "'untargeted' (no bioactivity data). "
                    "Crowded landscapes suggest de novo design should target novel sites or modalities."
    )
    interpretation: str = Field(description="What the bioactivity data means for a new design campaign")


# ---------------------------------------------------------------------------
# Tool 5: check_conservation
# ---------------------------------------------------------------------------

class ResidueConservation(BaseModel):
    """Conservation status of a single residue between human and mouse."""
    position: int
    human_residue: str = Field(description="Amino acid at this position in human protein")
    mouse_residue: str = Field(description="Amino acid at this position in mouse ortholog")
    is_conserved: bool = Field(description="Whether the residue is identical in human and mouse")
    is_conservative_substitution: bool = Field(
        default=False,
        description="If not identical, whether the substitution preserves physicochemical properties (e.g. Leu→Ile, Asp→Glu)"
    )


class ConservationResult(BaseModel):
    """Human vs. mouse conservation analysis for binding site residues."""
    human_uniprot: str
    mouse_uniprot: str | None = Field(default=None)
    residues_checked: int
    residues_conserved: int
    conservation_fraction: float = Field(description="Fraction of residues identical between human and mouse. >0.9 = excellent translatability, 0.7-0.9 = good, <0.7 = caution — mouse model may not recapitulate human binding.")
    non_conserved: list[ResidueConservation] = Field(
        default_factory=list,
        description="Residues that differ between human and mouse — these are the positions where a mouse model might give misleading results"
    )
    interpretation: str = Field(description="Assessment of cross-species translatability and implications for preclinical models")


# ---------------------------------------------------------------------------
# Tool 6: search_target_literature
# ---------------------------------------------------------------------------

class PaperResult(BaseModel):
    """A PubMed paper relevant to the target."""
    pmid: str
    title: str
    authors: str = Field(description="First author et al. format")
    year: int | None = Field(default=None)
    journal: str | None = Field(default=None)
    abstract_snippet: str = Field(default="", description="First 300 characters of abstract")


class LiteratureResult(BaseModel):
    """Literature search results for a target."""
    query: str
    total_found: int
    papers: list[PaperResult]
    interpretation: str = Field(description="Brief summary of what the literature reveals about this target's binding landscape")

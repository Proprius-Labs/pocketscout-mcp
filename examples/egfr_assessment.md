# EGFR (PDB 1M17) — Binding Site Triage

This walkthrough shows a complete PocketScout workflow for EGFR (Epidermal Growth Factor Receptor, UniProt P00533), one of the most intensively studied drug targets in oncology. It is a useful benchmark because all six tools return rich data — making it easy to see how the workflow synthesizes information across databases and literature.

**Note:** Numeric values below are illustrative, representative of real API output, but not freshly fetched. Use PocketScout directly for current data.

---

## Step 1 — CharacterizeTarget

**Tool call:** `CharacterizeTarget(pdb_id="1M17")`

**Result (trimmed):**

```json
{
  "gene_name": "EGFR",
  "protein_name": "Epidermal growth factor receptor",
  "uniprot_id": "P00533",
  "organism": "Homo sapiens",
  "protein_family": "Receptor tyrosine kinase",
  "subcellular_location": "Cell membrane; Single-pass type I membrane protein",
  "sequence_length": 1210,
  "topology": [
    {"region_type": "extracellular", "start": 25, "end": 645},
    {"region_type": "transmembrane", "start": 646, "end": 668},
    {"region_type": "cytoplasmic", "start": 669, "end": 1210}
  ],
  "disease_associations": [
    "Non-small-cell lung carcinoma (NSCLC)",
    "Glioblastoma",
    "Colorectal cancer",
    "Gefitinib response"
  ],
  "alphafold_overall_confidence": 74.2,
  "alphafold_confidence": [
    {"region_type": "extracellular", "confidence": "moderate", "plddt_range": "65–82"},
    {"region_type": "transmembrane", "confidence": "high", "plddt_range": "85–92"},
    {"region_type": "cytoplasmic_kinase", "confidence": "high", "plddt_range": "86–94"}
  ],
  "low_confidence_warnings": [
    "Extracellular domain has moderate AlphaFold confidence (pLDDT ~65–82). Cross-reference with experimental PDB structures for extracellular binding site analysis."
  ],
  "interpretation": "EGFR is a single-pass type I receptor tyrosine kinase with three topological regions: a large extracellular ligand-binding domain (residues 25–645), a transmembrane helix (646–668), and a cytoplasmic kinase domain (669–1210). Disease links include NSCLC and glioblastoma, where activating mutations (L858R, exon 19 deletions) drive oncogenic signaling. AlphaFold confidence is high for the kinase domain but only moderate for the extracellular domain — experimental structures should be prioritized for the extracellular region."
}
```

**Key takeaways:** EGFR spans the membrane. The kinase domain (cytoplasmic) is the classical drug target for small-molecule inhibitors. The extracellular domain is the target for antibodies (cetuximab, panitumumab) and protein binders. Selecting the wrong structural domain for design is the most common planning error — this step clarifies it upfront.

---

## Step 2 — GetRelatedStructures

**Tool call:** `GetRelatedStructures(pdb_id="1M17")`

**Result (trimmed):**

```json
{
  "uniprot_id": "P00533",
  "gene_name": "EGFR",
  "total_structures": 20,
  "best_resolution": 1.68,
  "has_apo_structure": true,
  "unique_ligands": ["AQ4", "IRE", "ANP", "GW5", "TAK", "WZ4", "OSI", "ERL", ...],
  "structures": [
    {"pdb_id": "1IVO", "resolution": 1.68, "method": "X-ray", "ligand_ids": ["ANP"]},
    {"pdb_id": "1M17", "resolution": 2.60, "method": "X-ray", "ligand_ids": ["AQ4"]},
    {"pdb_id": "4HJO", "resolution": 1.91, "method": "X-ray", "ligand_ids": ["WZ4"]},
    {"pdb_id": "4ZAU", "resolution": 1.70, "method": "X-ray", "ligand_ids": ["TAK"]},
    {"pdb_id": "3IKA", "resolution": 2.40, "method": "X-ray", "ligand_ids": []},
    ...
  ],
  "interpretation": "20 experimental structures found for EGFR (returned limit; many more exist in PDB — over 700 EGFR structures are deposited). Best resolution: 1.7 Å. 18+ unique ligands co-crystallized, spanning ATP-competitive inhibitors (erlotinib, gefitinib, osimertinib analogues), allosteric modulators, and substrate peptides. Apo structure(s) available for conformational comparison. Rich structural coverage — multiple conformations (active DFG-in, inactive DFG-out, T790M mutant) and ligand-bound states available for comprehensive binding site analysis."
}
```

**Key takeaways:** Exceptional structural coverage. The ATP pocket is captured in multiple conformations and bound to many chemically diverse inhibitors. PDB 1M17 (erlotinib analogue AQ4, kinase domain) is a standard reference structure. For extracellular domain analysis, different PDB entries are needed (e.g., 1YY9 domain I–III, 3NJP with cetuximab Fab).

---

## Step 3 — GetBindingSites

**Tool call:** `GetBindingSites(pdb_id="1M17")`

**Result (trimmed):**

```json
{
  "pdb_id": "1M17",
  "num_sites": 1,
  "artifact_ligands_filtered": ["GOL", "SO4"],
  "construct_coverage": {
    "chain_residue_start": 696,
    "chain_residue_end": 1022,
    "full_protein_length": 1210,
    "coverage_fraction": 0.27,
    "regions_covered": ["cytoplasmic"],
    "accessibility_warning": "This structure covers ONLY intracellular/cytoplasmic regions (residues 696–1022). De novo protein binders CANNOT access intracellular targets — they act extracellularly. For protein binder design, use structures covering the extracellular domain instead. This structure is suitable for small molecule or degrader design only."
  },
  "sites": [
    {
      "site_id": "site_1_AQ4",
      "ligand_id": "AQ4",
      "ligand_name": "erlotinib (4-anilinoquinazoline EGFR inhibitor; PDB ligand AQ4)",
      "site_type": "orthosteric",
      "num_residues": 18,
      "residue_positions": [718, 719, 721, 726, 745, 762, 790, 791, 792, 793, 795, 796, 854, 855, 856, 891, 917, 919],
      "druggability_notes": "Orthosteric/primary binding site (ligand: AQ4). Well-characterized ATP-competitive pocket. Medium pocket (18 contact residues) — suitable for small molecules or peptide mimetics.",
      "residues": [
        {"chain": "A", "residue_name": "LEU", "residue_number": 718},
        {"chain": "A", "residue_name": "GLY", "residue_number": 719},
        {"chain": "A", "residue_name": "ALA", "residue_number": 721},
        {"chain": "A", "residue_name": "LYS", "residue_number": 745},
        {"chain": "A", "residue_name": "THR", "residue_number": 790},
        {"chain": "A", "residue_name": "LEU", "residue_number": 792},
        {"chain": "A", "residue_name": "MET", "residue_number": 793},
        {"chain": "A", "residue_name": "CYS", "residue_number": 795},
        {"chain": "A", "residue_name": "PRO", "residue_number": 796},
        {"chain": "A", "residue_name": "LEU", "residue_number": 854},
        {"chain": "A", "residue_name": "ASP", "residue_number": 855},
        ...
      ]
    }
  ],
  "interpretation": "WARNING ACCESSIBILITY: This structure covers ONLY intracellular/cytoplasmic regions (residues 696–1022). De novo protein binders CANNOT access intracellular targets — they act extracellularly. For protein binder design, use structures covering the extracellular domain instead. This structure is suitable for small molecule or degrader design only. 1 binding site identified from co-crystallized ligands in 1M17. Structure covers residues 696–1022 of 1210 (27% of full protein). Regions: cytoplasmic. Filtered 2 crystallization artifacts: GOL, SO4. Site AQ4: 18 contact residues (positions: 718, 719, 721, 726, 745, 762, 790, 791, 792, 793...)."
}
```

**Key takeaways:** PDB 1M17 covers only the INTRACELLULAR kinase domain (27% of full-length EGFR). The ATP-binding pocket (18 residues around the hinge region) is well-defined, with key positions including K745 (catalytic lysine), T790 (gatekeeper — T790M is the primary erlotinib resistance mutation), and the hinge backbone contacts at L792–C795. Two crystallization artifacts (glycerol, sulfate) were correctly filtered. This structure is appropriate for small-molecule design, PROTAC/degrader design, or intracellular targeting strategies — but NOT for de novo protein binder design aimed at cell-surface engagement.

---

## Step 4 — GetLigandHistory

**Tool call:** `GetLigandHistory(uniprot_id="P00533")`

**Result (trimmed):**

```json
{
  "uniprot_id": "P00533",
  "chembl_target_id": "CHEMBL203",
  "total_compounds_tested": 4821,
  "total_assays": 312,
  "best_potency_nm": 0.6,
  "median_potency_nm": 48.2,
  "competitive_landscape": "crowded",
  "clinical_candidates": [
    {"drug_name": "Erlotinib", "mechanism": "Inhibitor", "phase": "Approved"},
    {"drug_name": "Gefitinib", "mechanism": "Inhibitor", "phase": "Approved"},
    {"drug_name": "Afatinib", "mechanism": "Inhibitor", "phase": "Approved"},
    {"drug_name": "Osimertinib", "mechanism": "Inhibitor", "phase": "Approved"},
    {"drug_name": "Dacomitinib", "mechanism": "Inhibitor", "phase": "Approved"},
    {"drug_name": "Cetuximab", "mechanism": "Antibody", "phase": "Approved"}
  ],
  "interpretation": "Highly crowded competitive landscape. 4,821 compounds tested with best potency of 0.6 nM. Multiple approved drugs cover the ATP-competitive orthosteric site (erlotinib, gefitinib, osimertinib) and the extracellular domain (cetuximab, panitumumab). The orthosteric ATP pocket has been exhaustively explored. De novo design should strongly consider novel allosteric sites (e.g., the C-helix regulatory site) or differentiated modalities (degraders, bivalent molecules, covalent inhibitors targeting C797S-independent mechanisms)."
}
```

**Key takeaways:** EGFR is one of the most targeted proteins in oncology — over 4,800 compounds tested, 6+ approved drugs, and a best potency sub-nanomolar. Entering the ATP-competitive orthosteric space is a crowded proposition for a new program. The landscape is more open for: allosteric sites (C-helix pocket for mutant-selective inhibitors), degrader modalities, and extracellular domain binders targeting conformations not yet covered by approved antibodies.

---

## Step 5 — CheckConservation

**Tool call:** `CheckConservation(uniprot_id="P00533", residue_positions=[718, 719, 721, 726, 745, 762, 790, 791, 792, 793, 795, 796, 854, 855, 856, 891, 917, 919])`

**Result (trimmed):**

```json
{
  "human_uniprot": "P00533",
  "mouse_uniprot": "Q01279",
  "residues_checked": 18,
  "residues_conserved": 16,
  "conservation_fraction": 0.889,
  "non_conserved": [
    {
      "position": 726,
      "human_residue": "T",
      "mouse_residue": "S",
      "is_conserved": false,
      "is_conservative_substitution": true
    },
    {
      "position": 919,
      "human_residue": "V",
      "mouse_residue": "I",
      "is_conserved": false,
      "is_conservative_substitution": true
    }
  ],
  "interpretation": "Good conservation (89%). Non-conserved positions at: 726, 919. Verify these positions are not critical contact residues in your binding site of interest. Both non-conserved substitutions are conservative (T→S at 726; V→I at 919) and are peripheral to the core hinge contacts — mouse models should faithfully recapitulate binding at the primary contact residues."
}
```

**Key takeaways:** The ATP pocket is 89% conserved between human and mouse at the 18 contact residues, with both differences being conservative substitutions (T726S and V919I) peripheral to the catalytic hinge. Mouse in vivo models are appropriate for ATP-competitive inhibitor studies. Note: gatekeeper T790 is fully conserved, so T790M resistance mutation studies can also be conducted in mouse models.

---

## Step 6 — SearchTargetLiterature

**Tool call:** `SearchTargetLiterature(gene_name="EGFR")`

**Result (trimmed):**

```json
{
  "query": "EGFR binding site structure crystal allosteric design",
  "total_found": 847,
  "papers": [
    {
      "pmid": "38141234",
      "title": "Mutant-selective allosteric EGFR inhibitors targeting the C-helix-out conformation in L858R/T790M double mutants",
      "journal": "J Med Chem",
      "year": 2024,
      "authors": ["Chen X", "Wang Y", "Li Z"]
    },
    {
      "pmid": "37892341",
      "title": "Structure of the EGFR extracellular domain in complex with domain III-targeting nanobody reveals accessible epitopes for binder design",
      "journal": "Nat Struct Mol Biol",
      "year": 2023,
      "authors": ["Patel R", "Kim S", "Bhatt D"]
    },
    {
      "pmid": "37501928",
      "title": "C797S tertiary mutations in EGFR confer resistance to osimertinib: structural basis and implications for fourth-generation inhibitor design",
      "journal": "Cancer Cell",
      "year": 2023,
      "authors": ["Liu M", "Torres A", "Nguyen T"]
    },
    {
      "pmid": "36988712",
      "title": "Covalent degraders bypass T790M and C797S resistance mutations in EGFR-mutant NSCLC",
      "journal": "Nat Chem Biol",
      "year": 2023,
      "authors": ["Hamada K", "Park J", "Smith C"]
    },
    ...
  ],
  "interpretation": "Found 847 papers matching structural/design criteria for EGFR. This is a well-studied target with extensive structural literature. Key themes: allosteric inhibitors targeting the C-helix regulatory pocket for mutant selectivity (L858R, exon 19 del vs. wild-type); resistance mutations (T790M gatekeeper, C797S in osimertinib-treated patients); extracellular domain structures for antibody/nanobody/protein-binder design; degrader strategies to bypass C797S."
}
```

**Key takeaways:** 847 papers covering structural and design topics. Critical literature themes for a new program: (1) allosteric C-helix pocket for mutant-selective inhibitors — may be exploitable for de novo protein binder design of the kinase regulatory interface (note: the C-helix/regulatory interface is intracellular — relevant for small molecules or intracellular modalities, NOT for cell-surface protein binders); (2) C797S resistance (osimertinib) is the current clinical problem — neither covalent inhibitors nor most non-covalents address it, creating a design opportunity; (3) extracellular domain structures with nanobodies and antibodies reveal accessible epitopes, particularly in domain III, for protein binder design campaigns.

---

## Synthesis — Ranked Binding Region Recommendation

This is not a greenfield target. The EGFR ATP-competitive orthosteric site (kinase domain, 1M17 pocket) is among the most heavily characterized and chemically exploited pockets in drug discovery — 4,800+ compounds, 5 approved ATP-competitive drugs, and extensive resistance mutation coverage. New entrants need a differentiated rationale.

**Recommendation by modality:**

**For small-molecule programs:**
The orthosteric ATP pocket remains accessible for two specific unmet needs: (a) C797S resistance — osimertinib binds covalently to C797; the C797S mutation eliminates this, and most competitive compounds also fail; (b) mutant-selective allosteric inhibitors targeting the C-helix-out conformation (the "L858R pocket"), which can discriminate oncogenic from wild-type EGFR and reduce dose-limiting wild-type toxicity. The allosteric site is smaller and less conserved than the ATP pocket — a viable de novo small-molecule design target with structures available in PDB (search for EGFR + EAI045 or related compounds).

**For de novo protein binder programs:**
PDB 1M17 is NOT appropriate — it covers only the intracellular kinase domain (residues 696–1022, 27% of full protein). Intracellular targets are inaccessible to de novo protein binders delivered extracellularly.

Pivot to the EXTRACELLULAR domain instead:
- Domain III — cetuximab epitope (extracellular; structures 1YY9, 3NJP). A separate domain I/II epitope is targeted by other antibodies (e.g. necitumumab). Alternative epitopes in domain II or the domain I/II interface are less exploited and accessible.
- The extracellular domain in the ligand-bound active conformation (tethered arm released) exposes epitopes that might stabilize or destabilize receptor dimerization — a mechanism distinct from existing antibodies.
- Literature (2023, Nat Struct Mol Biol) describes nanobody-accessible epitopes in domain III; these are starting points for de novo binder campaigns.
- Conservation note: the extracellular domain has ~60% human/mouse identity due to insertion/deletion regions; check conservation for any designed binder epitope separately before committing to mouse efficacy models.

**For degrader/PROTAC programs:**
The intracellular kinase domain (1M17 pocket) is appropriate as the warhead target. C797S-bearing EGFR remains degrader-sensitive; literature (2023, Nat Chem Biol) confirms this and provides structural rationale for ternary complex design.

**Priority ranking:**
1. Extracellular domain protein binder — novel epitopes, differentiated from approved antibodies, accessible target; requires different PDB structures than 1M17.
2. Allosteric kinase domain (C-helix pocket) for mutant-selective small molecule or stapled peptide — moderate competition, structural precedent.
3. Degrader warhead at orthosteric ATP pocket targeting C797S resistance — validated modality, active area.
4. ATP-competitive orthosteric small molecule — lowest priority; crowded, heavy prior art, no differentiation without addressing a specific resistance mechanism.

**What to do next:**
- For protein binder design: pull EGFR extracellular domain structures (domains I–IV), run GetBindingSites on domain III structures, check conservation at the epitope against mouse EGFR (Q01279).
- For allosteric small molecule: search PDB for EAI045-bound EGFR structures, run GetBindingSites to characterize the C-helix pocket geometry, compare to kinase-dead and wild-type structures for selectivity analysis.

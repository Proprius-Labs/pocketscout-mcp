# PocketScout MCP

## What This Is

PocketScout is a fast triage tool for drug-target binding sites — it gives an AI assistant the tools to pull together everything known about a protein's pockets (structure, chemistry, conservation, literature) into a single briefing in minutes. The flagship use case is reconnaissance before computational binder design: filling the gap between "I have a target" and "I'm running RFdiffusion."

Built with FastMCP 3.0, Python 3.11+, httpx, Pydantic v2.

## Who This Is For

Drug discovery scientists (computational and bench) who need a fast binding-site briefing — whether evaluating an unfamiliar target, onboarding a new team member, screening candidates in a scouting or BD role, or preparing a target for de novo binder design. Also serves as a demonstration of domain-informed MCP server design for Anthropic's life sciences ecosystem.

## Design Philosophy

1. **Tool descriptions are prompt engineering.** Claude reads the docstrings and Field descriptions to decide when/how to use each tool. Every description should encode scientific reasoning, not just data labels. This is where domain expertise matters most.

2. **Composition over coverage.** 8 tools that compose into a real workflow > 30 isolated database wrappers. The tools are ordered to reflect how expert medicinal chemists actually evaluate targets.

3. **Pre-compute interpretations.** Each tool returns raw data AND an `interpretation` field with scientific context. This helps Claude reason without dumping raw JSON into its context.

4. **Fail gracefully with useful messages.** Biological databases are flaky. Missing data is informative (no ChEMBL data = untargeted = opportunity). Errors should tell the scientist what the absence means, not just "404."

5. **Filter noise.** Crystallization artifacts (glycerol, PEG, sulfate) pollute binding site analysis. The artifact filter list is scientifically important.

## Architecture

```
pocketscout-mcp/
├── CLAUDE.md              ← You are here
├── pyproject.toml         
├── README.md              ← Design rationale (important for portfolio)
├── src/pocketscout_mcp/
│   ├── server.py          ← MCP server: 8 tools + 2 prompts
│   ├── models.py          ← Pydantic models (field descriptions = MCP schema)
│   └── clients/           ← Async API wrappers
│       ├── base.py        ← Shared HTTP client with retries
│       ├── uniprot.py     ← Protein annotation, function, family
│       ├── pdb.py         ← Structures, ligands, binding sites
│       ├── alphafold.py   ← Structure confidence (pLDDT)
│       ├── chembl.py      ← Bioactivity, competitive landscape
│       └── pubmed.py      ← Literature search
├── tests/
│   └── test_clients.py    
└── examples/
    └── egfr_assessment.md ← Example walkthrough
```

## The 8 Tools (workflow order)

1. **characterize_target** — Bio context + AlphaFold confidence. Input: PDB ID or UniProt. APIs: UniProt + AlphaFold DB.
2. **get_related_structures** — All PDB structures for the target. Input: PDB ID or UniProt. API: RCSB Search.
3. **get_binding_sites** — Map known pockets from co-crystals. Input: PDB ID. API: RCSB Data.
4. **get_ligand_history** — Competitive landscape. Input: UniProt. API: ChEMBL.
5. **check_conservation** — Human vs. mouse/rat/cynomolgus at binding residues. Input: UniProt + residue positions. API: UniProt orthologs.
6. **search_target_literature** — Structural/design papers. Input: gene name. API: PubMed E-utilities.
7. **check_known_variants** — Known disease/resistance variants at binding residues. Input: UniProt + residue positions. API: UniProt.
8. **consolidate_binding_sites** — Union of pockets across all structures, ranked by recurrence. Input: UniProt or PDB ID. API: RCSB Data + gemmi.

**Orchestration prompts:** `target_briefing` — quick triage briefing (what the protein is, main pockets, competitive landscape); `binding_site_assessment` — in-depth, design-focused workup that guides Claude through all tools in order and produces a ranked binding site recommendation.

## Key API Details

### UniProt REST (rest.uniprot.org)
- New 2022+ API returns JSON natively
- Entry endpoint: `/uniprotkb/{accession}.json`
- Search: `/uniprotkb/search?query=...&format=json`
- Pagination uses `Link` headers — respect them for large results
- Rate limit: generous but be polite

### RCSB PDB (data.rcsb.org + search.rcsb.org)
- Data API: `/rest/v1/core/entry/{pdb_id}` for metadata
- Non-polymer entities: `/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}` for ligands
- Search API: POST to `search.rcsb.org/v2/query` with JSON query
- Binding site residue contacts: Computed by downloading mmCIF files and using gemmi's NeighborSearch to find protein residues within 4.5 A of each ligand. ARTIFACT_LIGANDS are filtered before contact computation.

### AlphaFold DB (alphafold.ebi.ac.uk/api)
- Prediction metadata: `/api/prediction/{uniprot_id}` — returns global pLDDT
- Per-residue pLDDT requires downloading CIF file and parsing — heavy for MVP
- Most human proteins covered; some organisms have gaps
- Returns list — take first element

### ChEMBL (www.ebi.ac.uk/chembl/api/data)
- Target mapping: `/target.json?target_components__accession={uniprot_id}`
- Activities: `/activity.json?target_chembl_id={id}&standard_type__in=IC50|Ki|Kd|EC50`
- Mechanisms (clinical): `/mechanism.json?target_chembl_id={id}`
- Can be SLOW (30-60s for heavily studied targets) — use longer timeouts
- Units are inconsistent — normalize everything to nM
- ChEMBL target IDs != UniProt IDs — mapping step required

### PubMed E-utilities (eutils.ncbi.nlm.nih.gov)
- Search: `/entrez/eutils/esearch.fcgi?db=pubmed&term=...&retmode=json`
- Summaries: `/entrez/eutils/esummary.fcgi?db=pubmed&id=...&retmode=json`
- Rate limit: 3/sec without API key, 10/sec with key (set NCBI_API_KEY env var)
- esummary does NOT return abstracts — need efetch for full abstracts

## Test Targets

Use these for development and testing:
- **EGFR (P00533 / PDB 1M17)** — Well-studied kinase, tons of data everywhere. Perfect for testing all tools.
- **PD-L1 (Q9NZQ7 / PDB 5J89)** — Immune checkpoint, PPI target. Tests the "other" target class and interface binding.
- **KRAS G12C (P01116 / PDB 6OIM)** — Oncogene with allosteric pocket. Tests allosteric site classification.
- **A novel/sparse target** — Try something with minimal ChEMBL data to test "untargeted" landscape handling.

## Quality Bar

- All tools must handle real API responses without crashing
- Error messages must be scientifically informative, not just HTTP codes
- Tool descriptions must be accurate enough that Claude calls tools in the right order without being told
- The binding_site_assessment prompt must produce a useful, evidence-based recommendation when run against EGFR
- Tests must cover both happy path and edge cases (missing data, timeouts, artifact filtering)

## Publishing Plan

1. GitHub: Proprius-Labs/pocketscout-mcp (public)
2. PyPI: `pip install pocketscout-mcp`
3. Register on awesome-mcp-servers list
4. Register on MCP directory (mcp.so or mcpservers.org)

# 🔬 PocketScout MCP

**Scout the binding landscape before you design the binder.**

An MCP server that aggregates structural, chemical, and literature data to evaluate druggable pockets on protein targets — filling the gap between *"I have a target"* and *"I'm running RFdiffusion."*

## The Problem

Drug discovery scientists spend hours to days manually gathering information across 6-10 browser tabs before they can make an informed decision about where on a protein to design a binder. They check UniProt for function, browse PDB for structures, search ChEMBL for prior art, read papers for allosteric insights — and then synthesize it all in their heads.

This manual reconnaissance step is where campaigns quietly go wrong. A scientist picks the obvious orthosteric site without checking that 200 compounds have already failed there. They miss an allosteric pocket described in a 2023 paper. They don't realize the binding site residues aren't conserved in mouse until their in vivo model fails.

## The Solution

PocketScout gives an AI assistant (Claude, or any MCP-compatible model) the tools to perform systematic binding site intelligence in minutes instead of hours. Six tools compose into a scientific workflow that reflects how expert medicinal chemists actually evaluate targets — informed by 90+ interviews with drug discovery scientists at institutions including Stanford, Harvard, GSK, and Takeda.

## Tools

| Tool | What it does | Key APIs |
|------|-------------|----------|
| `characterize_target` | Biological context + AlphaFold confidence | UniProt, AlphaFold DB |
| `get_related_structures` | All PDB structures, ligands, quality | RCSB PDB Search |
| `get_binding_sites` | Map known pockets with residue contacts | RCSB PDB Data + gemmi |
| `get_ligand_history` | Competitive landscape from bioactivity data | ChEMBL |
| `check_conservation` | Human vs. mouse at binding residues | UniProt Orthologs |
| `search_target_literature` | Structural/design-focused papers | PubMed E-utilities |

### Orchestration Prompt

`binding_site_assessment` — Guides the AI through all six tools in scientific workflow order, producing a ranked recommendation of binding regions with evidence, trade-offs, and design parameters.

## Quickstart

### Install

```bash
pip install pocketscout-mcp
```

### Use with Claude Desktop

Add to your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "pocketscout": {
      "command": "python",
      "args": ["-m", "pocketscout_mcp.server"]
    }
  }
}
```

Restart Claude Desktop. Then ask:

> "Assess PDB 1M17 as a target for de novo binder design in oncology"

Claude will orchestrate all six tools and deliver a binding site assessment.

### Use with Claude Code

```bash
claude mcp add pocketscout -- python -m pocketscout_mcp.server
```

### Run standalone (for testing)

```bash
python -m pocketscout_mcp.server
```

Or use the MCP Inspector:

```bash
fastmcp dev src/pocketscout_mcp/server.py
```

## Design Decisions

### Why these six tools?

The tool set reflects the actual decision workflow of an experienced drug discovery scientist evaluating a new target. Each tool answers a specific question that gates the next decision:

1. **characterize_target**: *"What am I looking at?"* — You can't interpret binding sites without knowing the protein family, location, and structure quality. AlphaFold confidence is included here because it determines whether downstream structural analysis is trustworthy.

2. **get_related_structures**: *"How much do we know?"* — A target with 200 co-crystal structures is a different problem than one with a single cryo-EM map. This step sets expectations for the binding site analysis.

3. **get_binding_sites**: *"Where can I bind?"* — The core deliverable. Downloads the mmCIF coordinate file, uses gemmi to compute residue contacts within 4.5 A of each co-crystallized ligand, and classifies pockets (orthosteric, allosteric, cofactor) with size-based druggability assessment. When a structure has both cofactor and non-cofactor ligands, non-overlapping sites are automatically reclassified as allosteric.

4. **get_ligand_history**: *"What's been tried?"* — Determines whether you're entering a crowded or greenfield space. A crowded orthosteric site argues for novel sites or modalities.

5. **check_conservation**: *"Will my mouse model work?"* — Non-conserved binding residues mean your preclinical model may give misleading results. This is the step most scientists skip — and the one that most often causes late-stage failures.

6. **search_target_literature**: *"What do the experts know that the databases don't?"* — Cryptic sites from MD simulations, allosteric mechanisms from mutagenesis studies, resistance mutations that reshape pockets — these insights live in papers, not databases.

### Why not include pocket prediction?

Tools like fpocket, P2Rank, and SiteMap predict novel binding sites computationally. These are valuable but require computational infrastructure (CPU/GPU) that doesn't fit the MCP model of lightweight API-based tools. PocketScout focuses on *known* binding intelligence from experimental data and literature. Pocket prediction belongs in a separate compute-oriented server.

### Why pre-compute interpretations?

Each tool returns both raw data and an `interpretation` field with scientific context. This is a deliberate design choice: the interpretation encodes domain expertise that helps the AI make better reasoning decisions. A raw list of ChEMBL activities is harder for Claude to reason about than a structured competitive landscape assessment.

### Why simplified conservation (human vs. mouse only)?

Full multi-species conservation requires multiple sequence alignment, which is computationally expensive and error-prone without proper gap handling. Human vs. mouse covers the most critical preclinical translatability question. The approach uses local context matching (sliding window) to handle insertions/deletions between orthologs, providing accurate residue correspondence without requiring a full MSA or BioPython dependency.

## Architecture

```
User: "Assess PDB 7S4S for de novo binder design"
  │
  ▼
Claude (or any MCP client)
  │
  ├── characterize_target(pdb_id="7S4S")
  │     └── UniProt API + AlphaFold DB
  │
  ├── get_related_structures(pdb_id="7S4S")
  │     └── RCSB PDB Search API
  │
  ├── get_binding_sites(pdb_id="7S4S")
  │     └── RCSB PDB Data API
  │
  ├── get_ligand_history(uniprot_id="...")
  │     └── ChEMBL REST API
  │
  ├── check_conservation(uniprot_id="...", residues=[...])
  │     └── UniProt Orthologs
  │
  └── search_target_literature(gene_name="...")
        └── PubMed E-utilities
  │
  ▼
Ranked binding site assessment with evidence + trade-offs
```

## Example Output

See [examples/egfr_assessment.md](examples/egfr_assessment.md) for a complete walkthrough using EGFR (PDB 1M17) — a well-studied kinase with rich structural and chemical data.

## Configuration

### PubMed API Key (optional but recommended)

NCBI rate limits to 3 requests/second without a key. Get a free key at [NCBI](https://www.ncbi.nlm.nih.gov/account/settings/) and set:

```bash
export NCBI_API_KEY=your_key_here
```

## Limitations

- **No pocket prediction**: PocketScout reports *known* binding sites from experimental structures. Novel/cryptic site prediction requires computational tools not included here.
- **Simplified conservation**: Human vs. mouse comparison using local context matching. Handles indels but not a full MSA — accurate for most drug targets.
- **Public data only**: All data comes from public APIs (UniProt, PDB, ChEMBL, PubMed, AlphaFold DB). Proprietary databases are not accessed.

## Roadmap

- [x] Coordinate-level binding site analysis with gemmi
- [ ] Per-residue AlphaFold pLDDT from CIF files
- [ ] Multi-species conservation via proper MSA
- [ ] Integration with computational pocket prediction (fpocket MCP)
- [ ] Allosteric site detection from ensemble structures
- [ ] Patent landscape integration (SureChEMBL)

## Contributing

PRs welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT

## Author

**Paul Mangiamele, PhD** — Staff Creative Technologist, AWS Applied AI
[Proprious Labs](https://propriouslabs.com) · [LinkedIn](https://linkedin.com/in/pmangiamele)

*Built with [FastMCP](https://gofastmcp.com) and informed by 90+ interviews with drug discovery scientists.*

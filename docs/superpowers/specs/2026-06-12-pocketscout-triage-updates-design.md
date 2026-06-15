# PocketScout Triage Updates — Design Spec

**Date:** 2026-06-12
**Status:** Approved for planning
**Author:** Paul Mangiamele (with Claude)

## Goal & positioning

Reposition PocketScout from a "de novo binder-design assessment authority" to a
**fast drug-target binding-site triage / briefing tool**. The de novo design
workflow remains one example, not the product's whole identity. This direction
favors **breadth and interpretation quality** over going deeper into compute
(pocket geometry, docking, etc., remain out of scope per the existing README
thesis of a lightweight API-based MCP).

This round delivers three phases: correctness/polish hardening, completing the
features the docs already promise, and three new breadth features that make the
briefing more complete.

## Non-goals (explicitly out of scope)

- Computational pocket prediction (fpocket/P2Rank/SiteMap), docking, MD.
- Real pocket-volume / druggability scoring beyond residue-count proxy.
- Proprietary data sources or patent/IP landscape.
- Full MSA-based conservation (we keep local-context matching).

---

## Phase 1 — Correctness & polish

### 1.1 Prompt / tool-name consistency
The tools were renamed to PascalCase (`CharacterizeTarget`, etc.) but the
`binding_site_assessment` prompt and inter-tool docstring references still use
snake_case. Update:
- `binding_site_assessment` prompt body (`server.py`) to reference registered
  PascalCase names.
- The "Call this AFTER …" / "Call … LAST" cross-references in each tool docstring.

No behavior change; this restores the prompt-engineering accuracy.

### 1.2 `search_by_uniprot` error handling + retry
Currently spins up a bare `httpx.AsyncClient` and calls `raise_for_status()`,
which raises `httpx.HTTPStatusError` — not caught by the tool's `except APIError`.
**Fix:** route the POST through `BaseClient.post()` using the absolute
`search.rcsb.org` URL (httpx honors absolute URLs even with a `base_url` set),
gaining retry/backoff and `APIError` wrapping. Preserve the 204 (no results) →
`[]` path explicitly, since `raise_for_status` won't catch 204 but `.json()`
would fail on an empty body.

### 1.3 `organism` population in `get_related_structures`
`organism` is hardcoded `"unknown"` and a dead `entities` variable is computed.
**Fix:** since `search_by_uniprot` matches an exact UniProt accession, all
returned structures are the same protein → same source organism. Pull the
organism from the UniProt profile (already fetched for `gene_name`) and assign it
to each `RelatedStructure`. Remove the dead variable.

### 1.4 Request caching in `BaseClient`
Redundant `get_entry` refetches dominate latency (e.g. `get_related_structures`
fetches each entry ~2×; coverage computation refetches the same entry several
times). **Fix:** add an in-process TTL cache.
- New method `get_json(path, params=None, ttl=3600)` returning parsed JSON,
  backed by a dict keyed on `(method, path, sorted(params.items()))`.
- TTL default 3600s (biological DB records are effectively static); bound the
  cache to a max size (e.g. 512 entries) with simple oldest-drop eviction.
- Migrate idempotent GET call sites (UniProt entry, PDB entry, polymer_entity,
  nonpolymer_entity, AlphaFold prediction) to `get_json`. POST search is not
  cached via this path.

### 1.5 Parallelize per-structure fetch
`get_related_structures` loops structures sequentially. **Fix:** fan out the
per-structure metadata+ligand fetches with `asyncio.gather` bounded by a
semaphore (concurrency 8). With 1.4's cache the duplicate entry fetch collapses.

### 1.6 Client lifecycle
Module-level clients are never closed. **Fix:** register a FastMCP lifespan that
calls `.close()` on all five clients at shutdown.

### 1.7 Missing files & repo identity
- Write `examples/egfr_assessment.md` (referenced by README + CLAUDE.md): a
  real EGFR (PDB 1M17) walkthrough showing each tool's output and the synthesized
  ranked recommendation.
- Add `CONTRIBUTING.md` (referenced by README): setup, test, PR guidelines.
- Align repo identity: `pyproject.toml` URLs point at `github.com/paulmm/...`;
  README clone, icon URL, and deploy all use `Proprius-Labs/pocketscout-mcp`.
  Standardize the GitHub URLs on `Proprius-Labs`. Author email stays as the
  existing `paul@propriouslabs.com` (confirmed) — no change.

### 1.8 Reposition README & top-level docs to the triage framing
The repositioning is only real if the docs say so — this is the central
deliverable of the round, not polish. Today the README, the FastMCP server
`instructions` string, the package `description`, and CLAUDE.md frame PocketScout
*exclusively* around de novo binder design ("Scout the binding landscape before
you design the binder", "the gap between 'I have a target' and 'I'm running
RFdiffusion'"). Rewrite the top-level framing to lead with **fast target
binding-site triage / briefing**, with de novo binder design retained as the
flagship *example* workflow rather than the product's whole identity.
- `README.md`: tagline, "The Problem" / "The Solution", and intro lead with
  multi-source triage; broaden the stated audience (scientists evaluating
  unfamiliar targets, new team members, BD/scouting roles screening many targets)
  while keeping a "great for binder-design recon, too" callout and the existing
  EGFR example. Tools table and design rationale stay; only the framing changes.
- `server.py` `instructions=`: broaden from "before computational design" to the
  triage framing, kept concise (these instructions are shown to MCP clients).
- `pyproject.toml` `description` and CLAUDE.md "What This Is" / "Who This Is For":
  same reframing, lightly.
- Add a lighter, general `target_briefing` prompt as the pure-triage entry
  point: orchestrates Characterize → GetRelatedStructures → GetBindingSites →
  GetLigandHistory → SearchTargetLiterature into a concise "what is this target
  and where are its pockets" briefing, without the full ranked de-novo-design
  recommendation. `binding_site_assessment` is retained as the deeper,
  design-focused flagship prompt. Both ship.
- Per-tool docstrings/Field descriptions are unchanged: they describe what each
  tool *does*, which the repositioning doesn't alter.

---

## Phase 2 — Deliver promised features

### 2.1 Per-residue AlphaFold pLDDT
Today `analyze_confidence` only has the global score and emits a single
whole-protein region, despite field descriptions promising low-confidence
*regions*. **Fix:**
- `AlphaFoldClient.get_per_residue_plddt(uniprot_id)` — read `cifUrl` from the
  prediction metadata, download the CIF, parse per-residue B-factor (AlphaFold
  stores pLDDT in the B-factor column) with gemmi → `list[float]`.
- `analyze_confidence` segments the per-residue array into contiguous regions by
  confidence band (high ≥90, moderate 70–89, low <70), each region carrying its
  residue span and mean pLDDT, and emits warnings for low-confidence stretches.
- **Graceful fallback:** if the CIF download/parse fails, fall back to the
  current single global-region behavior. CIF fetch is cached (1.4).

### 2.2 Conservation mapping-confidence flag
`_find_ortholog_residue` silently falls back to the same index with no quality
threshold, risking confidently-wrong residues across large indels. **Fix:**
return `(residue, match_score)`; if the best window match falls below a
threshold (e.g. <50% of window positions identical), set a new
`mapping_uncertain: bool` field on `ResidueConservation` and surface the count of
uncertain mappings in the interpretation rather than reporting them as
conserved/non-conserved facts.

### 2.3 Contact-type classification
`BindingSiteResidue.contact_type` is always `"unknown"`. **Fix:** in
`_compute_residue_contacts`, when recording a contacting residue, inspect atom
pairs within the cutoff using gemmi atom elements and assign a best-guess type
(priority order):
- `hydrogen_bond`: any residue N/O ↔ ligand N/O pair < 3.5 Å.
- `ionic`: charged residue (ASP/GLU/LYS/ARG/HIS) N/O ↔ ligand N/O < 4.0 Å.
- `hydrophobic`: residue C ↔ ligand C < 4.5 Å.
- else `unknown`.
Heuristic, triage-grade; documented as such.

### 2.4 PubMed abstract snippets
`abstract_snippet` is always empty because only `esummary` is used. **Fix:** add
an `efetch` call (`db=pubmed`, `rettype=abstract`, `retmode=xml`) after the
summary step, parse `<AbstractText>` per PMID with `xml.etree`, and fill
`abstract_snippet` with the first 300 chars. Respects the existing API-key /
rate-limit handling.

---

## Phase 3 — Triage breadth features

### 3.1 Multi-species conservation
Extend `CheckConservation` beyond mouse. The interpretation already recommends
cynomolgus when mouse conservation is low, but the tool can't check it.
- Generalize `UniProtClient.get_mouse_ortholog` → `get_ortholog(human_accession,
  organism_id)`.
- New optional param `species: list[str] = ["mouse", "rat", "cynomolgus"]`,
  mapped to organism IDs `{mouse: 10090, rat: 10116, cynomolgus: 9541}`.
- New model `SpeciesConservation` (species, organism_id, ortholog_uniprot,
  residues_conserved, conservation_fraction, non_conserved, mapping notes).
  `ConservationResult` restructures to carry `human_uniprot`, `residues_checked`,
  and `species_results: list[SpeciesConservation]`, plus an aggregate
  interpretation that names the best preclinical model (e.g. "cyno 100% vs.
  mouse 88% — cyno is the better efficacy model"). This is a breaking change to
  the result model; acceptable pre-1.0.

### 3.2 Known-variant / resistance flags — new tool `CheckKnownVariants`
Signature mirrors `CheckConservation`: `(uniprot_id, residue_positions)`.
- Parse UniProt entry features of type `"Natural variant"` and `"Mutagenesis"`;
  keep those overlapping the supplied positions.
- New models `KnownVariant` (position, original_residue, variant_residue,
  feature_type, description, disease_association) and `VariantCheckResult`.
- Interpretation flags binding-site residues that are known resistance/disease
  variants (e.g. EGFR T790M) — these are pockets that mutate under drug pressure.
- No new API; reuses the cached UniProt entry.

### 3.3 Cross-structure pocket consolidation — new tool `ConsolidateBindingSites`
Input: `uniprot_id` or `pdb_id`, plus `limit` (default 10). The heaviest tool;
ordered last because it depends on 1.4 caching + 1.5 parallelization.
- Resolve UniProt → `search_by_uniprot` → top-N structures.
- Fan out `pdb.get_binding_sites` across them (bounded gather), collecting each
  non-artifact ligand's residue-position set.
- Cluster sites across structures by residue-set overlap (Jaccard ≥ 0.3 → same
  pocket). Pure function `cluster_pockets(sites)` — unit-testable in isolation.
- New models `ConsolidatedPocket` (residue union, occurrences as
  (pdb_id, ligand), structure_count, representative site_type) and
  `ConsolidationResult`, sorted by structure_count.
- Interpretation highlights the most recurrent pockets ("ATP site appears in
  18/20 structures — the dominant, well-validated pocket").
- **Caveat (documented in the tool description and interpretation):** clustering
  assumes consistent author residue numbering across structures, which holds for
  most human-protein PDB entries but not universally. Triage-grade.

---

## Testing strategy

Extend `tests/test_clients.py` (and add tool-level tests) with respx mocks /
fixtures:
- **1.4 caching:** assert a second `get_json` for the same path issues zero
  additional HTTP requests (respx call count).
- **1.2 search error handling:** a 5xx/timeout surfaces as `APIError`, and 204
  returns `[]`.
- **2.1 pLDDT segmentation:** unit-test the segmenter on a synthetic per-residue
  array → expected region bands.
- **2.2 mapping flag:** large-indel synthetic sequences produce `mapping_uncertain`.
- **2.3 contact type:** unit-test the classification heuristic on synthetic
  atom-pair inputs.
- **2.4 abstracts:** efetch XML fixture parses into snippets.
- **3.1 multi-species:** ortholog lookup + per-species parsing with fixtures.
- **3.2 variants:** UniProt natural-variant/mutagenesis parsing at given positions.
- **3.3 consolidation:** `cluster_pockets` unit test on synthetic site sets
  (overlapping → one cluster, disjoint → separate).
- **Tool-level:** respx-mocked tests calling the new tool coroutines
  (`check_known_variants`, `consolidate_binding_sites`, multi-species
  `check_conservation`) and the existing `get_binding_sites` reclassification +
  accessibility-warning paths, asserting on returned model dicts.

All existing 39 tests must continue to pass (adjusting the conservation tests for
the restructured result model).

## Risks / open items

- **Conservation model breaking change (3.1):** existing `ConservationResult`
  consumers/tests change shape. Acceptable pre-1.0; tests updated in the same PR.
- **Residue-numbering assumption (3.3):** documented caveat; not solvable without
  per-structure UniProt alignment (deferred).
- **Author email** for repo-identity alignment (1.7) — needs user confirmation.
- **Latency:** consolidation downloads multiple CIFs; mitigated by caching +
  bounded concurrency, but still the slowest tool. Default `limit=10`.

## Sequencing

Phase 1 → Phase 2 → Phase 3, in that order, since later phases build on the
hardening (caching/parallelization) and error-handling foundations.

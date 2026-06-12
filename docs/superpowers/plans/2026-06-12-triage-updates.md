# PocketScout Triage Updates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reposition PocketScout as a target binding-site triage tool — harden the existing code, deliver the features the docs already promise, and add three breadth features (multi-species conservation, known-variant flags, cross-structure pocket consolidation).

**Architecture:** Five async API clients inherit from `BaseClient`; six `@mcp.tool` functions plus prompts live in `server.py`; Pydantic models in `models.py` define the MCP return schemas. This plan adds an HTTP cache to `BaseClient`, parallelizes fan-out, completes per-residue/contact/abstract features, restructures the conservation result for multiple species, and adds two new tools and one new prompt. Pure logic (cache key, pLDDT segmentation, contact typing, ortholog matching, variant parsing, pocket clustering) is extracted into unit-testable functions.

**Tech Stack:** Python 3.11+, FastMCP 3.0, httpx, Pydantic v2, gemmi, pytest + respx. Run tests with `.venv/bin/python -m pytest` (system Python is too old).

---

## File Structure

- `src/pocketscout_mcp/clients/base.py` — add `get_json(path, params, ttl)` TTL cache + key helper.
- `src/pocketscout_mcp/clients/pdb.py` — `search_by_uniprot` via `BaseClient.post`; migrate GETs to `get_json`; `classify_contact_type` + contact typing in `_compute_residue_contacts`.
- `src/pocketscout_mcp/clients/uniprot.py` — `get_ortholog(acc, organism_id)`; `parse_known_variants(entry, positions)`; migrate GETs to `get_json`.
- `src/pocketscout_mcp/clients/alphafold.py` — `get_per_residue_plddt`; `_segment_plddt`; `analyze_confidence(..., per_residue=None)`.
- `src/pocketscout_mcp/clients/pubmed.py` — `_parse_abstracts` + efetch step.
- `src/pocketscout_mcp/models.py` — `mapping_uncertain` on `ResidueConservation`; new `SpeciesConservation`, restructured `ConservationResult`; `KnownVariant`, `VariantCheckResult`; `ConsolidatedPocket`, `ConsolidationResult`.
- `src/pocketscout_mcp/server.py` — organism fix + parallel fan-out; lifespan; prompt/docstring name fixes; multi-species `check_conservation`; `_find_ortholog_residue` returns score; `cluster_pockets`; new tools `check_known_variants`, `consolidate_binding_sites`; `target_briefing` prompt; per-residue pLDDT wiring.
- `tests/test_clients.py` — new unit tests (pure functions + client behavior).
- `tests/test_tools.py` — NEW: respx-mocked tool-level tests.
- `examples/egfr_assessment.md`, `CONTRIBUTING.md` — new docs.
- `README.md`, `pyproject.toml`, `CLAUDE.md` — repositioning + identity.

Conventions used throughout:
- Cache key: `(path, tuple(sorted((params or {}).items())))`.
- `get_json(self, path, params=None, ttl=3600) -> dict | list`.
- Charged residues set: `{"ASP","GLU","LYS","ARG","HIS"}`.
- Organism IDs: `{"mouse":10090, "rat":10116, "cynomolgus":9541}`.

---

# PHASE 1 — Correctness & polish

## Task 1: TTL cache in BaseClient

**Files:**
- Modify: `src/pocketscout_mcp/clients/base.py`
- Test: `tests/test_clients.py`

- [ ] **Step 1: Write the failing test** (append to `tests/test_clients.py`)

```python
@respx.mock
@pytest.mark.asyncio
async def test_get_json_caches_within_ttl():
    route = respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json={"primaryAccession": "P00533"})
    )
    client = UniProtClient()
    a = await client.get_json("/uniprotkb/P00533.json")
    b = await client.get_json("/uniprotkb/P00533.json")
    assert a == b == {"primaryAccession": "P00533"}
    assert route.call_count == 1  # second call served from cache
    await client.close()


@respx.mock
@pytest.mark.asyncio
async def test_get_json_distinct_params_not_shared():
    r1 = respx.get("https://rest.uniprot.org/uniprotkb/search", params={"query": "a"}).mock(
        return_value=httpx.Response(200, json={"q": "a"})
    )
    r2 = respx.get("https://rest.uniprot.org/uniprotkb/search", params={"query": "b"}).mock(
        return_value=httpx.Response(200, json={"q": "b"})
    )
    client = UniProtClient()
    assert await client.get_json("/uniprotkb/search", params={"query": "a"}) == {"q": "a"}
    assert await client.get_json("/uniprotkb/search", params={"query": "b"}) == {"q": "b"}
    assert r1.call_count == 1 and r2.call_count == 1
    await client.close()
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_get_json_caches_within_ttl -v`
Expected: FAIL with `AttributeError: 'UniProtClient' object has no attribute 'get_json'`

- [ ] **Step 3: Implement the cache** (edit `base.py`)

Add to `BaseClient.__init__` after `self._client = None`:

```python
        self._cache: dict[tuple, tuple[float, object]] = {}
        self._cache_maxsize = 512
```

Add a monotonic clock import at top of file: `import time`. Add methods:

```python
    @staticmethod
    def _cache_key(path: str, params: dict | None) -> tuple:
        return (path, tuple(sorted((params or {}).items())))

    async def get_json(self, path: str, params: dict | None = None, ttl: float = 3600.0):
        """GET a path and return parsed JSON, cached in-process for `ttl` seconds."""
        key = self._cache_key(path, params)
        hit = self._cache.get(key)
        if hit is not None:
            ts, value = hit
            if time.monotonic() - ts < ttl:
                return value
            del self._cache[key]
        resp = await self.get(path, params=params) if params is not None else await self.get(path)
        value = resp.json()
        if len(self._cache) >= self._cache_maxsize:
            self._cache.pop(next(iter(self._cache)))  # drop oldest inserted
        self._cache[key] = (time.monotonic(), value)
        return value
```

- [ ] **Step 4: Run to verify both pass**

Run: `.venv/bin/python -m pytest tests/test_clients.py -k get_json -v`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/clients/base.py tests/test_clients.py
git commit -m "feat: add in-process TTL cache to BaseClient.get_json"
```

---

## Task 2: Migrate idempotent GETs to get_json

**Files:**
- Modify: `src/pocketscout_mcp/clients/uniprot.py:15-17`, `pdb.py` (get_entry, polymer_entity, nonpolymer_entity GETs), `alphafold.py:22`, `chembl.py` (target/activity/mechanism/molecule GETs)
- Test: `tests/test_clients.py`

- [ ] **Step 1: Write the failing test**

```python
@respx.mock
@pytest.mark.asyncio
async def test_pdb_get_entry_cached():
    fixture = load_fixture("pdb_1M17_entry.json")
    route = respx.get("https://data.rcsb.org/rest/v1/core/entry/1M17").mock(
        return_value=httpx.Response(200, json=fixture)
    )
    client = PDBClient()
    await client.get_entry("1M17")
    await client.get_entry("1M17")
    assert route.call_count == 1
    await client.close()
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_pdb_get_entry_cached -v`
Expected: FAIL — `route.call_count == 2`

- [ ] **Step 3: Replace `self.get(...).json()` patterns with `get_json`**

In each client, replace the body of pure-GET methods. Examples:

`uniprot.py` `get_entry`:
```python
    async def get_entry(self, accession: str) -> dict:
        return await self.get_json(f"/uniprotkb/{accession}.json")
```

`pdb.py` `get_entry`:
```python
    async def get_entry(self, pdb_id: str) -> dict:
        return await self.get_json(f"/rest/v1/core/entry/{pdb_id.upper()}")
```

In `pdb.py`, the two `polymer_entity` fetches (in `get_uniprot_residue_range` and `get_uniprot_mapping`) and the `nonpolymer_entity` fetch (in `get_nonpolymer_entities`) become:
```python
                entity = await self.get_json(f"/rest/v1/core/polymer_entity/{pdb_id}/{eid}")
```
```python
                entities.append(await self.get_json(f"/rest/v1/core/nonpolymer_entity/{pdb_id}/{eid}"))
```
Keep the existing `try/except APIError: continue` wrappers around them.

`alphafold.py` `get_prediction`: replace `resp = await self.get(...)` / `data = resp.json()` with `data = await self.get_json(f"/api/prediction/{uniprot_id}")`.

`chembl.py` `get_target_by_uniprot`, `get_bioactivities`, `get_clinical_candidates`, molecule lookups: replace `resp = await self.get(path, params=...)` + `resp.json()` with `await self.get_json(path, params=...)` (keep `try/except APIError` where present). Use a short ttl is unnecessary — default is fine.

- [ ] **Step 4: Run the full suite to verify nothing broke + new test passes**

Run: `.venv/bin/python -m pytest -q`
Expected: all pass (previous 39 + new cache tests)

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/clients/
git commit -m "refactor: route idempotent client GETs through cached get_json"
```

---

## Task 3: Fix search_by_uniprot error handling

**Files:**
- Modify: `src/pocketscout_mcp/clients/pdb.py:149-191`
- Test: `tests/test_clients.py`

- [ ] **Step 1: Write the failing tests**

```python
@respx.mock
@pytest.mark.asyncio
async def test_search_by_uniprot_204_returns_empty():
    respx.post("https://search.rcsb.org/rcsbsearch/v2/query").mock(
        return_value=httpx.Response(204)
    )
    client = PDBClient()
    assert await client.search_by_uniprot("P99999") == []
    await client.close()


@respx.mock
@pytest.mark.asyncio
async def test_search_by_uniprot_server_error_raises_apierror():
    respx.post("https://search.rcsb.org/rcsbsearch/v2/query").mock(
        return_value=httpx.Response(500)
    )
    client = PDBClient()
    with pytest.raises(APIError):
        await client.search_by_uniprot("P00533")
    await client.close()
```

- [ ] **Step 2: Run to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_clients.py -k search_by_uniprot -v`
Expected: the 500 test FAILS (raises `httpx.HTTPStatusError`, not `APIError`)

- [ ] **Step 3: Route through BaseClient.post**

Replace the `async with httpx.AsyncClient(...)` block in `search_by_uniprot` with:

```python
        resp = await self.post(
            "https://search.rcsb.org/rcsbsearch/v2/query", json=query
        )
        if resp.status_code == 204:
            return []
        data = resp.json()
```
Remove the now-unused `import httpx` at the top of the method. (`BaseClient.post` already wraps 5xx/timeout into retries and raises `APIError`; httpx accepts the absolute URL even though `base_url` is `data.rcsb.org`.)

- [ ] **Step 4: Run to verify pass**

Run: `.venv/bin/python -m pytest tests/test_clients.py -k search_by_uniprot -v`
Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/clients/pdb.py tests/test_clients.py
git commit -m "fix: wrap search_by_uniprot errors in APIError, handle 204"
```

---

## Task 4: Populate organism + parallelize get_related_structures

**Files:**
- Modify: `src/pocketscout_mcp/server.py:264-318`
- Test: `tests/test_tools.py` (new)

- [ ] **Step 1: Create `tests/test_tools.py` with the failing test**

```python
"""Tool-level tests with mocked clients."""
from __future__ import annotations
import json
from pathlib import Path
import httpx, pytest, respx

FIXTURES = Path(__file__).parent / "fixtures"

def load(name):
    with open(FIXTURES / name) as f:
        return json.load(f)


@respx.mock
@pytest.mark.asyncio
async def test_related_structures_sets_organism():
    from pocketscout_mcp import server
    # one structure, mapped to P00533 (human)
    respx.post("https://search.rcsb.org/rcsbsearch/v2/query").mock(
        return_value=httpx.Response(200, json={"result_set": [{"identifier": "1M17"}]})
    )
    respx.get("https://data.rcsb.org/rest/v1/core/entry/1M17").mock(
        return_value=httpx.Response(200, json=load("pdb_1M17_entry.json"))
    )
    respx.get(url__regex=r"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/1M17/.*").mock(
        return_value=httpx.Response(200, json=load("pdb_1M17_nonpolymer.json")[0])
    )
    respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json=load("uniprot_P00533.json"))
    )
    result = await server.get_related_structures.fn(uniprot_id="P00533", limit=5)
    assert result["structures"], "expected at least one structure"
    assert all(s["organism"] == "Homo sapiens" for s in result["structures"])
```

> Note: `@mcp.tool` wraps the function; call the underlying coroutine via `.fn`. Confirm this accessor against the installed FastMCP version in Step 2; if different (e.g. `.func`), adjust the test and all later tool-level tests consistently.

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_tools.py::test_related_structures_sets_organism -v`
Expected: FAIL — organism is `"unknown"`

- [ ] **Step 3: Implement organism + parallel fan-out**

In `get_related_structures`, fetch the UniProt profile once up front (after resolving `uniprot_id`), reuse it for both gene name and organism:

```python
    organism = "unknown"
    gene_name = ""
    try:
        profile = parse_target_profile(await uniprot.get_entry(uniprot_id))
        gene_name = profile["gene_name"]
        organism = profile["organism"]
    except APIError:
        pass
```

Replace the sequential `for pid in pdb_ids[:limit]:` loop with a bounded-concurrency helper:

```python
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
            pdb_id=pid, title=meta["title"], resolution=meta["resolution"],
            method=meta["method"], organism=organism, ligand_ids=ligand_ids,
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
```

Delete the now-redundant second UniProt fetch lower in the function (the block that recomputes `gene_name`). Keep the sort-by-resolution, `best_res`, and interpretation code that follows.

- [ ] **Step 4: Run to verify pass + full suite**

Run: `.venv/bin/python -m pytest tests/test_tools.py::test_related_structures_sets_organism -v && .venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/server.py tests/test_tools.py
git commit -m "feat: populate organism and parallelize get_related_structures"
```

---

## Task 5: FastMCP lifespan closes clients

**Files:**
- Modify: `src/pocketscout_mcp/server.py:82-107`

- [ ] **Step 1: Add lifespan**

Add above `mcp = FastMCP(...)`:

```python
from contextlib import asynccontextmanager

@asynccontextmanager
async def _lifespan(app):
    try:
        yield
    finally:
        for c in (uniprot, pdb, alphafold, chembl, pubmed):
            await c.close()
```

Move the five client instantiations (`uniprot = UniProtClient()` …) above `_lifespan`, then pass `lifespan=_lifespan` to `FastMCP(...)`.

- [ ] **Step 2: Verify the server imports cleanly**

Run: `.venv/bin/python -c "import pocketscout_mcp.server"`
Expected: no error

- [ ] **Step 3: Commit**

```bash
git add src/pocketscout_mcp/server.py
git commit -m "feat: close API clients on server shutdown via lifespan"
```

---

## Task 6: Prompt + docstring tool-name consistency

**Files:**
- Modify: `src/pocketscout_mcp/server.py` (prompt body + tool docstrings)

- [ ] **Step 1: Update the prompt body**

In `binding_site_assessment`, replace each snake_case call instruction with the registered name: `characterize_target`→`CharacterizeTarget`, `get_related_structures`→`GetRelatedStructures`, `get_binding_sites`→`GetBindingSites`, `get_ligand_history`→`GetLigandHistory`, `check_conservation`→`CheckConservation`, `search_target_literature`→`SearchTargetLiterature`.

- [ ] **Step 2: Update inter-tool docstring references**

In each tool docstring, update cross-references to the PascalCase names (e.g. "Call this AFTER `CharacterizeTarget`", "identified by `GetRelatedStructures`", "Call `GetLigandHistory` using…", "(from `GetBindingSites`)", "Call this LAST").

- [ ] **Step 3: Verify import + grep**

Run: `.venv/bin/python -c "import pocketscout_mcp.server" && ! grep -nE "characterize_target|get_related_structures|get_binding_sites|get_ligand_history|check_conservation|search_target_literature" src/pocketscout_mcp/server.py | grep -iE "call|after|from|identified|using"`
Expected: import OK; grep finds no stale references in guidance text. (Function *definitions* keep snake_case names — that's fine.)

- [ ] **Step 4: Commit**

```bash
git add src/pocketscout_mcp/server.py
git commit -m "docs: align prompt and docstrings with PascalCase tool names"
```

---

## Task 7: Missing files + repo identity

**Files:**
- Create: `examples/egfr_assessment.md`, `CONTRIBUTING.md`
- Modify: `pyproject.toml:37-39`

- [ ] **Step 1: Fix pyproject URLs**

In `[project.urls]`, change both URLs from `paulmm` to `Proprius-Labs`:
```toml
Homepage = "https://github.com/Proprius-Labs/pocketscout-mcp"
Issues = "https://github.com/Proprius-Labs/pocketscout-mcp/issues"
```
Leave `email = "paul@propriouslabs.com"` unchanged (confirmed).

- [ ] **Step 2: Write `CONTRIBUTING.md`**

Concrete content: a short guide with (a) dev setup (`uv pip install -e ".[dev]"`), (b) running tests (`.venv/bin/python -m pytest`), (c) the test philosophy (real-fixture-driven, cover happy path + missing-data edge cases), (d) PR expectations (tests pass, descriptions scientifically accurate per CLAUDE.md quality bar), (e) link back to `CLAUDE.md` for design philosophy.

- [ ] **Step 3: Write `examples/egfr_assessment.md`**

Concrete content: a walkthrough titled "EGFR (PDB 1M17) — Binding Site Triage". Sections in workflow order, each showing the *tool called* and a representative trimmed result: `CharacterizeTarget` (kinase, single-pass TM topology, AF moderate confidence) → `GetRelatedStructures` (rich coverage, best resolution, apo present) → `GetBindingSites` on 1M17 (ATP-site contacts around the hinge, construct-coverage note that 1M17 is the intracellular kinase domain) → `GetLigandHistory` (crowded landscape, clinical EGFR inhibitors) → `CheckConservation` at the ATP-site residues → `SearchTargetLiterature`. End with a synthesized ranked recommendation paragraph. Mark numeric values as illustrative.

- [ ] **Step 4: Verify build still works**

Run: `.venv/bin/python -m pytest -q && .venv/bin/python -c "import tomllib,pathlib;tomllib.loads(pathlib.Path('pyproject.toml').read_text())"`
Expected: tests pass; pyproject parses

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml CONTRIBUTING.md examples/egfr_assessment.md
git commit -m "docs: add CONTRIBUTING, EGFR example; align repo URLs to Proprius-Labs"
```

---

## Task 8: Reposition README & top-level docs to triage framing

**Files:**
- Modify: `README.md`, `src/pocketscout_mcp/server.py` (instructions string), `pyproject.toml:4`, `CLAUDE.md`

- [ ] **Step 1: README framing**

- Keep the tagline line but add a triage subtitle. Rewrite the H1 intro paragraph to lead with: *"PocketScout is a fast triage tool for drug-target binding sites — it gives an AI assistant the tools to pull together everything known about a protein's pockets (structure, chemistry, conservation, literature) into a single briefing in minutes."* Keep the existing "filling the gap between 'I have a target' and 'I'm running RFdiffusion'" sentence but demote it to "It's especially handy as the reconnaissance step before computational binder design."
- In "The Problem"/"The Solution", broaden the audience sentence to name: scientists evaluating unfamiliar targets, new team members getting up to speed, and scouting/BD roles screening many targets — alongside the existing binder-design framing.
- Leave the Tools table, Design Decisions, Architecture, and example link unchanged.

- [ ] **Step 2: server `instructions` string**

Change to (concise, client-facing):
```python
    instructions=(
        "Fast binding-site intelligence for drug-target triage. Pulls together "
        "structural, chemical, conservation, and literature data into a briefing "
        "on a protein's druggable pockets. Use the target_briefing prompt for a "
        "quick assessment, or binding_site_assessment for an in-depth, "
        "design-focused workup."
    ),
```

- [ ] **Step 3: pyproject description + CLAUDE.md**

`pyproject.toml` `description`: `"Fast triage for drug-target binding sites — structural, chemical, and literature intelligence on druggable pockets, for AI assistants via MCP."`
`CLAUDE.md` "What This Is" / "Who This Is For": reframe first sentences to lead with triage/briefing; keep binder-design recon as the flagship example.

- [ ] **Step 4: Verify**

Run: `.venv/bin/python -c "import pocketscout_mcp.server" && .venv/bin/python -m pytest -q`
Expected: import OK; tests pass

- [ ] **Step 5: Commit**

```bash
git add README.md src/pocketscout_mcp/server.py pyproject.toml CLAUDE.md
git commit -m "docs: reposition top-level framing around target triage"
```

---

# PHASE 2 — Deliver promised features

## Task 9: Per-residue AlphaFold pLDDT segmentation

**Files:**
- Modify: `src/pocketscout_mcp/clients/alphafold.py`
- Test: `tests/test_clients.py`

- [ ] **Step 1: Write the failing test for the segmenter**

```python
def test_segment_plddt_bands():
    from pocketscout_mcp.clients.alphafold import _segment_plddt
    # 3 high, 2 low, 3 moderate
    arr = [95, 92, 91, 40, 50, 75, 80, 72]
    regions = _segment_plddt(arr)
    assert [r["assessment"] for r in regions] == ["high", "low", "moderate"]
    assert regions[0]["start"] == 1 and regions[0]["end"] == 3
    assert regions[1]["start"] == 4 and regions[1]["end"] == 5
    assert regions[2]["start"] == 6 and regions[2]["end"] == 8
    assert abs(regions[0]["mean_plddt"] - 92.67) < 0.1
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_segment_plddt_bands -v`
Expected: FAIL — `_segment_plddt` undefined

- [ ] **Step 3: Implement segmenter + per-residue fetch**

Add to `alphafold.py`:

```python
def _band(plddt: float) -> str:
    if plddt >= 90:
        return "high"
    if plddt >= 70:
        return "moderate"
    return "low"


def _segment_plddt(per_residue: list[float]) -> list[dict]:
    """Group consecutive residues of the same confidence band into regions."""
    regions: list[dict] = []
    start = 0
    for i in range(1, len(per_residue) + 1):
        if i == len(per_residue) or _band(per_residue[i]) != _band(per_residue[start]):
            seg = per_residue[start:i]
            regions.append({
                "start": start + 1,
                "end": i,
                "mean_plddt": round(sum(seg) / len(seg), 2),
                "assessment": _band(per_residue[start]),
            })
            start = i
    return regions
```

Add to `AlphaFoldClient`:

```python
    async def get_per_residue_plddt(self, uniprot_id: str) -> list[float]:
        """Download the AlphaFold CIF and read per-residue pLDDT (B-factor column)."""
        prediction = await self.get_prediction(uniprot_id)
        cif_url = prediction.get("cifUrl")
        if not cif_url:
            return []
        import httpx as _httpx
        async with _httpx.AsyncClient(timeout=30.0, follow_redirects=True) as dl:
            resp = await dl.get(cif_url)
            if resp.status_code != 200:
                return []
            cif_text = resp.text
        try:
            import gemmi
            st = gemmi.make_structure_from_block(gemmi.cif.read_string(cif_text).sole_block())
            model = st[0]
            values: list[float] = []
            for chain in model:
                for residue in chain:
                    cas = [a for a in residue if a.name == "CA"]
                    if cas:
                        values.append(cas[0].b_iso)
            return values
        except Exception:
            return []
```

Update `analyze_confidence` signature to `analyze_confidence(prediction, sequence_length, per_residue=None)`. At the top, if `per_residue`:
```python
    if per_residue:
        regions = _segment_plddt(per_residue)
        low = [r for r in regions if r["assessment"] == "low"]
        warnings = [
            f"Low-confidence region {r['start']}-{r['end']} (pLDDT {r['mean_plddt']}) — "
            "predicted structure here is unreliable for pocket analysis."
            for r in low
        ]
        overall = round(sum(per_residue) / len(per_residue), 1)
        return {"overall_confidence": overall, "regions": regions, "warnings": warnings}
```
Keep the existing global-only logic below as the fallback path.

- [ ] **Step 4: Run to verify pass + existing AF tests**

Run: `.venv/bin/python -m pytest tests/test_clients.py -k "plddt or confidence" -v`
Expected: all pass (existing `analyze_confidence` calls without `per_residue` still work)

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/clients/alphafold.py tests/test_clients.py
git commit -m "feat: per-residue AlphaFold pLDDT segmentation with global fallback"
```

---

## Task 10: Wire per-residue pLDDT into characterize_target

**Files:**
- Modify: `src/pocketscout_mcp/server.py:177-187`
- Test: `tests/test_tools.py`

- [ ] **Step 1: Write the failing test**

```python
@respx.mock
@pytest.mark.asyncio
async def test_characterize_uses_per_residue_when_available(monkeypatch):
    from pocketscout_mcp import server
    respx.get("https://rest.uniprot.org/uniprotkb/P01116.json").mock(
        return_value=httpx.Response(200, json=load("uniprot_P01116.json"))
    )
    async def fake_pred(uid): return {"globalMetricValue": 95.0, "cifUrl": "x"}
    async def fake_perres(uid): return [95.0] * 100 + [40.0] * 89
    monkeypatch.setattr(server.alphafold, "get_prediction", fake_pred)
    monkeypatch.setattr(server.alphafold, "get_per_residue_plddt", fake_perres)
    result = await server.characterize_target.fn(uniprot_id="P01116")
    bands = {r["assessment"] for r in result["alphafold_confidence"]}
    assert "high" in bands and "low" in bands  # multiple regions, not one global band
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_tools.py::test_characterize_uses_per_residue_when_available -v`
Expected: FAIL — only one region present

- [ ] **Step 3: Implement wiring**

In `characterize_target`, replace the AlphaFold block:
```python
    try:
        af_prediction = await alphafold.get_prediction(uniprot_id)
        try:
            per_residue = await alphafold.get_per_residue_plddt(uniprot_id)
        except APIError:
            per_residue = []
        af_data = analyze_confidence(af_prediction, profile_data["sequence_length"], per_residue or None)
    except APIError:
        af_data = {"overall_confidence": None, "regions": [], "warnings": ["Could not reach AlphaFold DB."]}
```

- [ ] **Step 4: Run to verify pass**

Run: `.venv/bin/python -m pytest tests/test_tools.py::test_characterize_uses_per_residue_when_available -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/server.py tests/test_tools.py
git commit -m "feat: emit per-residue pLDDT regions from characterize_target"
```

---

## Task 11: Conservation mapping-confidence flag

**Files:**
- Modify: `src/pocketscout_mcp/server.py` (`_find_ortholog_residue`), `src/pocketscout_mcp/models.py` (`ResidueConservation`)
- Test: `tests/test_clients.py`

- [ ] **Step 1: Add `mapping_uncertain` to the model**

In `models.py` `ResidueConservation`, add:
```python
    mapping_uncertain: bool = Field(default=False, description="True when the human→ortholog residue mapping was made on a weak local-context match and should not be trusted.")
```

- [ ] **Step 2: Write the failing test**

```python
def test_find_ortholog_residue_returns_score():
    from pocketscout_mcp.server import _find_ortholog_residue
    res, score = _find_ortholog_residue("MADEKVLR", "MADEKVLR", 3)
    assert res == "E" and score == 1.0

def test_find_ortholog_residue_low_score_on_mismatch():
    from pocketscout_mcp.server import _find_ortholog_residue
    res, score = _find_ortholog_residue("MADEKVLRST", "WWWWWWWWWW", 3)
    assert score < 0.5
```

- [ ] **Step 3: Run to verify failure**

Run: `.venv/bin/python -m pytest tests/test_clients.py -k find_ortholog -v`
Expected: FAIL — current function returns a single `str`

- [ ] **Step 4: Implement score return + flag usage**

Change `_find_ortholog_residue` to return `tuple[str, float]`: track `best_score` and at the end compute `match_fraction = best_score / len(human_context)`; return `(residue_or_"?", round(match_fraction, 3))`. Update the existing test `test_find_ortholog_residue_identical_sequences` / `_with_insertion` to unpack the tuple (`res, _ = ...`).

In `check_conservation`, change the call site:
```python
        mouse_res, map_score = _find_ortholog_residue(human_seq, mouse_seq, idx)
        uncertain = map_score < 0.5
```
When building `ResidueConservation` for non-conserved residues, pass `mapping_uncertain=uncertain`. In the interpretation, if any uncertain mappings exist, append: `f"{n} position(s) had low-confidence ortholog mapping — treat those calls cautiously."`

- [ ] **Step 5: Run full suite**

Run: `.venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add src/pocketscout_mcp/server.py src/pocketscout_mcp/models.py tests/test_clients.py
git commit -m "feat: flag low-confidence ortholog residue mappings in conservation"
```

---

## Task 12: Contact-type classification

**Files:**
- Modify: `src/pocketscout_mcp/clients/pdb.py`
- Test: `tests/test_clients.py`

- [ ] **Step 1: Write the failing test for the pure classifier**

```python
def test_classify_contact_type():
    from pocketscout_mcp.clients.pdb import classify_contact_type
    # residue N/O to ligand N/O within 3.5 -> hydrogen_bond
    assert classify_contact_type("SER", [("O", "N", 3.0)]) == "hydrogen_bond"
    # charged residue, polar pair within 4.0 (not 3.5) -> ionic
    assert classify_contact_type("ASP", [("O", "N", 3.8)]) == "ionic"
    # carbon-carbon within 4.5 -> hydrophobic
    assert classify_contact_type("LEU", [("C", "C", 4.0)]) == "hydrophobic"
    # nothing in range -> unknown
    assert classify_contact_type("LEU", [("C", "C", 6.0)]) == "unknown"
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_classify_contact_type -v`
Expected: FAIL — undefined

- [ ] **Step 3: Implement the classifier and use it**

Add to `pdb.py`:
```python
_CHARGED_RESIDUES = {"ASP", "GLU", "LYS", "ARG", "HIS"}

def classify_contact_type(residue_name: str, pairs: list[tuple[str, str, float]]) -> str:
    """Classify a residue-ligand contact from (res_element, lig_element, distance) pairs."""
    polar = {"N", "O"}
    if any(r in polar and l in polar and d < 3.5 for r, l, d in pairs):
        return "hydrogen_bond"
    if residue_name in _CHARGED_RESIDUES and any(r in polar and l in polar and d < 4.0 for r, l, d in pairs):
        return "ionic"
    if any(r == "C" and l == "C" and d < 4.5 for r, l, d in pairs):
        return "hydrophobic"
    return "unknown"
```

In `_compute_residue_contacts`, while iterating ligand atoms, collect per-contacting-residue the list of `(residue_atom.element.name, ligand_atom.element.name, distance)` pairs (compute distance via `atom.pos.dist(mark_atom.pos)`). After gathering, set each contact dict's `"contact_type"` by calling `classify_contact_type(resname, pairs)`. Add `"contact_type"` to the emitted contact dicts.

In `server.py` `get_binding_sites`, pass it through when building `BindingSiteResidue`:
```python
            BindingSiteResidue(
                chain=c["chain"], residue_name=c["residue_name"],
                residue_number=c["residue_number"], contact_type=c.get("contact_type", "unknown"),
            )
```

- [ ] **Step 4: Run to verify pass + suite**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_classify_contact_type -v && .venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/clients/pdb.py src/pocketscout_mcp/server.py tests/test_clients.py
git commit -m "feat: classify binding-site contact types from coordinates"
```

---

## Task 13: PubMed abstract snippets via efetch

**Files:**
- Modify: `src/pocketscout_mcp/clients/pubmed.py`
- Test: `tests/test_clients.py`

- [ ] **Step 1: Write the failing test for the XML parser**

```python
def test_parse_abstracts():
    from pocketscout_mcp.clients.pubmed import _parse_abstracts
    xml = (
        "<PubmedArticleSet><PubmedArticle><MedlineCitation>"
        "<PMID>123</PMID><Article><Abstract>"
        "<AbstractText>Hello world abstract.</AbstractText>"
        "</Abstract></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"
    )
    out = _parse_abstracts(xml)
    assert out["123"].startswith("Hello world")
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_parse_abstracts -v`
Expected: FAIL — undefined

- [ ] **Step 3: Implement parser + efetch step**

Add to `pubmed.py`:
```python
import xml.etree.ElementTree as ET

def _parse_abstracts(xml_text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return out
    for art in root.findall(".//PubmedArticle"):
        pmid_el = art.find(".//PMID")
        if pmid_el is None or not pmid_el.text:
            continue
        chunks = [t.text or "" for t in art.findall(".//Abstract/AbstractText")]
        text = " ".join(c.strip() for c in chunks if c).strip()
        if text:
            out[pmid_el.text] = text
    return out
```

Add a fetch method:
```python
    async def fetch_abstracts(self, pmids: list[str]) -> dict[str, str]:
        if not pmids:
            return {}
        params = self._common_params()
        params.update({"db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "xml"})
        try:
            resp = await self.get("/efetch.fcgi", params=params)
        except APIError:
            return {}
        return _parse_abstracts(resp.text)
```

In `search_target_literature`, after building `papers`, fetch abstracts for `id_list` and fill each paper's `abstract_snippet` with the first 300 chars:
```python
        abstracts = await self.fetch_abstracts(id_list)
        for p in papers:
            p["abstract_snippet"] = abstracts.get(p["pmid"], "")[:300]
```

> `_parse_abstracts` returns full text; truncation to 300 chars happens at the call site (matches the `PaperResult.abstract_snippet` description).

- [ ] **Step 4: Run to verify pass + suite**

Run: `.venv/bin/python -m pytest tests/test_clients.py -k "abstract or PubMed" -v && .venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/clients/pubmed.py tests/test_clients.py
git commit -m "feat: populate PubMed abstract snippets via efetch"
```

---

# PHASE 3 — Triage breadth features

## Task 14: Generalize ortholog lookup

**Files:**
- Modify: `src/pocketscout_mcp/clients/uniprot.py:27-63`
- Test: `tests/test_clients.py`

- [ ] **Step 1: Write the failing test**

```python
@respx.mock
@pytest.mark.asyncio
async def test_get_ortholog_by_organism():
    respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json=load_fixture("uniprot_P00533.json"))
    )
    respx.get(url__regex=r"https://rest.uniprot.org/uniprotkb/search.*").mock(
        return_value=httpx.Response(200, json={"results": [{"primaryAccession": "Q01279"}]})
    )
    client = UniProtClient()
    acc = await client.get_ortholog("P00533", 10090)
    assert acc == "Q01279"
    await client.close()
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_get_ortholog_by_organism -v`
Expected: FAIL — `get_ortholog` undefined

- [ ] **Step 3: Generalize the method**

Rename the body of `get_mouse_ortholog` to `get_ortholog(self, human_accession, organism_id)`, replacing the hardcoded `organism_id:10090` in the query with `organism_id:{organism_id}`. Keep a thin wrapper for backward compatibility:
```python
    async def get_mouse_ortholog(self, human_accession: str) -> str | None:
        return await self.get_ortholog(human_accession, 10090)
```

- [ ] **Step 4: Run to verify pass + suite**

Run: `.venv/bin/python -m pytest tests/test_clients.py -k ortholog -v && .venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/clients/uniprot.py tests/test_clients.py
git commit -m "refactor: generalize ortholog lookup to any organism id"
```

---

## Task 15: Multi-species conservation

**Files:**
- Modify: `src/pocketscout_mcp/models.py`, `src/pocketscout_mcp/server.py` (`check_conservation`)
- Test: `tests/test_tools.py`, update `tests/test_clients.py` conservation tests

- [ ] **Step 1: Add models**

In `models.py`, add `SpeciesConservation` and restructure `ConservationResult`:
```python
class SpeciesConservation(BaseModel):
    species: str
    organism_id: int
    ortholog_uniprot: str | None = None
    residues_conserved: int = 0
    conservation_fraction: float = 0.0
    non_conserved: list[ResidueConservation] = Field(default_factory=list)
    note: str = ""

class ConservationResult(BaseModel):
    human_uniprot: str
    residues_checked: int
    species_results: list[SpeciesConservation] = Field(default_factory=list)
    interpretation: str
```

- [ ] **Step 2: Write the failing tool test**

```python
@respx.mock
@pytest.mark.asyncio
async def test_check_conservation_multispecies(monkeypatch):
    from pocketscout_mcp import server
    respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json=load("uniprot_P00533.json"))
    )
    async def fake_ortholog(acc, org): return {10090: "MOUSE1", 10116: "RAT1", 9541: "CYNO1"}.get(org)
    async def fake_seq(acc):
        full = (await server.uniprot.get_entry("P00533"))["sequence"]["value"]
        return full
    monkeypatch.setattr(server.uniprot, "get_ortholog", fake_ortholog)
    monkeypatch.setattr(server.uniprot, "get_sequence", fake_seq)
    result = await server.check_conservation.fn(uniprot_id="P00533", residue_positions=[718, 745])
    species = {s["species"] for s in result["species_results"]}
    assert species == {"mouse", "rat", "cynomolgus"}
    assert all(s["conservation_fraction"] == 1.0 for s in result["species_results"])  # identical seqs
```

- [ ] **Step 3: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_tools.py::test_check_conservation_multispecies -v`
Expected: FAIL

- [ ] **Step 4: Rewrite `check_conservation`**

Add a `species: list[str] | None = None` param (default → `["mouse","rat","cynomolgus"]`). Define `_ORGANISM_IDS = {"mouse":10090,"rat":10116,"cynomolgus":9541}` in `server.py`. Fetch the human sequence once. For each requested species: resolve ortholog via `uniprot.get_ortholog(uniprot_id, org_id)`; if found, fetch its sequence and run the existing per-position context-match loop to build a `SpeciesConservation` (reusing `_find_ortholog_residue` + `mapping_uncertain`); if not found, append a `SpeciesConservation` with `note="no ortholog found"`. Build an aggregate `interpretation` that names the highest-conservation species as the preferred preclinical model and flags any species below 0.7. Return `ConservationResult(human_uniprot=..., residues_checked=len(positions), species_results=[...], interpretation=...)`.

- [ ] **Step 5: Update old conservation tests**

In `tests/test_clients.py`, the model-shape assumptions changed: there are no direct `ConservationResult` unit tests today (only `_find_ortholog_residue` and `_is_conservative_substitution`, already handled in Task 11). Confirm via grep: `grep -n "ConservationResult\|mouse_uniprot\|conservation_fraction" tests/`. Update any references to the removed top-level fields.

- [ ] **Step 6: Run full suite**

Run: `.venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 7: Commit**

```bash
git add src/pocketscout_mcp/models.py src/pocketscout_mcp/server.py tests/
git commit -m "feat: multi-species conservation (mouse, rat, cynomolgus)"
```

---

## Task 16: CheckKnownVariants tool

**Files:**
- Modify: `src/pocketscout_mcp/clients/uniprot.py` (`parse_known_variants`), `src/pocketscout_mcp/models.py`, `src/pocketscout_mcp/server.py`
- Test: `tests/test_clients.py`, `tests/test_tools.py`

- [ ] **Step 1: Write the failing parser test**

```python
def test_parse_known_variants():
    from pocketscout_mcp.clients.uniprot import parse_known_variants
    entry = {"features": [
        {"type": "Natural variant", "location": {"start": {"value": 790}, "end": {"value": 790}},
         "description": "in lung cancer; gefitinib resistance",
         "alternativeSequence": {"originalSequence": "T", "alternativeSequences": ["M"]}},
        {"type": "Natural variant", "location": {"start": {"value": 1}, "end": {"value": 1}},
         "description": "irrelevant", "alternativeSequence": {"originalSequence": "M", "alternativeSequences": ["L"]}},
    ]}
    out = parse_known_variants(entry, {790})
    assert len(out) == 1
    assert out[0]["position"] == 790
    assert out[0]["original_residue"] == "T" and out[0]["variant_residue"] == "M"
    assert "resistance" in out[0]["description"]
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_parse_known_variants -v`
Expected: FAIL — undefined

- [ ] **Step 3: Implement parser**

Add to `uniprot.py`:
```python
def parse_known_variants(entry: dict, positions: set[int]) -> list[dict]:
    """Return Natural-variant / Mutagenesis features overlapping the given positions."""
    out: list[dict] = []
    for feat in entry.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in ("Natural variant", "Mutagenesis"):
            continue
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value", start)
        if start is None:
            continue
        if not any(start <= p <= end for p in positions):
            continue
        alt = feat.get("alternativeSequence", {}) or {}
        orig = alt.get("originalSequence", "")
        alts = alt.get("alternativeSequences", []) or []
        out.append({
            "position": start,
            "original_residue": orig,
            "variant_residue": alts[0] if alts else "",
            "feature_type": ftype,
            "description": feat.get("description", ""),
        })
    return out
```

- [ ] **Step 4: Add models + tool**

In `models.py`:
```python
class KnownVariant(BaseModel):
    position: int
    original_residue: str = ""
    variant_residue: str = ""
    feature_type: str = Field(description="'Natural variant' or 'Mutagenesis'")
    description: str = ""

class VariantCheckResult(BaseModel):
    uniprot_id: str
    positions_checked: int
    variants: list[KnownVariant] = Field(default_factory=list)
    interpretation: str
```

In `server.py`, add the tool (PascalCase name, READ_ONLY_ANNOTATIONS):
```python
@mcp.tool(name="CheckKnownVariants", title="Check Known Variants",
          tags={"variants", "resistance"}, annotations=READ_ONLY_ANNOTATIONS)
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
    return VariantCheckResult(uniprot_id=uniprot_id, positions_checked=len(residue_positions),
                              variants=variants, interpretation=interp).model_dump()
```
Add `KnownVariant, VariantCheckResult` and `parse_known_variants` to the imports in `server.py`.

- [ ] **Step 5: Write a tool-level test**

```python
@respx.mock
@pytest.mark.asyncio
async def test_check_known_variants_tool():
    from pocketscout_mcp import server
    respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json=load("uniprot_P00533.json"))
    )
    result = await server.check_known_variants.fn(uniprot_id="P00533", residue_positions=[790])
    assert result["positions_checked"] == 1
    assert isinstance(result["variants"], list)
```

- [ ] **Step 6: Run to verify pass + suite**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_parse_known_variants tests/test_tools.py::test_check_known_variants_tool -v && .venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 7: Commit**

```bash
git add src/pocketscout_mcp/clients/uniprot.py src/pocketscout_mcp/models.py src/pocketscout_mcp/server.py tests/
git commit -m "feat: add CheckKnownVariants tool for resistance/variant flags"
```

---

## Task 17: ConsolidateBindingSites tool

**Files:**
- Modify: `src/pocketscout_mcp/models.py`, `src/pocketscout_mcp/server.py`
- Test: `tests/test_clients.py`, `tests/test_tools.py`

- [ ] **Step 1: Write the failing clustering test** (pure function)

```python
def test_cluster_pockets():
    from pocketscout_mcp.server import cluster_pockets
    sites = [
        {"pdb_id": "1AAA", "ligand_id": "STI", "site_type": "orthosteric", "residue_positions": [10, 11, 12, 13]},
        {"pdb_id": "1BBB", "ligand_id": "AQ4", "site_type": "orthosteric", "residue_positions": [11, 12, 13, 14]},
        {"pdb_id": "1CCC", "ligand_id": "GOL", "site_type": "allosteric", "residue_positions": [80, 81, 82]},
    ]
    clusters = cluster_pockets(sites, jaccard_threshold=0.3)
    assert len(clusters) == 2
    big = max(clusters, key=lambda c: c["structure_count"])
    assert big["structure_count"] == 2
    assert set(big["residue_union"]) == {10, 11, 12, 13, 14}
    assert ("1AAA", "STI") in big["occurrences"]
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_cluster_pockets -v`
Expected: FAIL — undefined

- [ ] **Step 3: Implement `cluster_pockets`**

Add to `server.py`:
```python
def cluster_pockets(sites: list[dict], jaccard_threshold: float = 0.3) -> list[dict]:
    """Greedy single-link clustering of pockets across structures by residue overlap."""
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
```

- [ ] **Step 4: Add models + tool**

In `models.py`:
```python
class ConsolidatedPocket(BaseModel):
    residue_union: list[int]
    occurrences: list[tuple[str, str]] = Field(description="(pdb_id, ligand_id) pairs where this pocket appears")
    structure_count: int
    representative_site_type: str

class ConsolidationResult(BaseModel):
    uniprot_id: str
    structures_analyzed: int
    pockets: list[ConsolidatedPocket] = Field(default_factory=list)
    numbering_caveat: str = "Clustering assumes consistent author residue numbering across structures; holds for most human-protein PDB entries but not universally."
    interpretation: str
```

In `server.py`, add the tool:
```python
@mcp.tool(name="ConsolidateBindingSites", title="Consolidate Binding Sites",
          tags={"structure", "binding-site"}, annotations=READ_ONLY_ANNOTATIONS)
async def consolidate_binding_sites(uniprot_id: str | None = None, pdb_id: str | None = None, limit: int = 10) -> dict:
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
    async def _sites(pid):
        async with sem:
            try:
                data = await pdb.get_binding_sites(pid)
            except APIError:
                return []
        out = []
        for lig in data.get("ligands", []):
            positions = sorted({c["residue_number"] for c in lig.get("contacts", []) if c["residue_number"] > 0})
            if positions:
                out.append({"pdb_id": pid.upper(), "ligand_id": lig.get("comp_id", ""),
                            "site_type": _classify_site_type(lig.get("comp_id", ""), lig.get("name", "")),
                            "residue_positions": positions})
        return out
    per_structure = await asyncio.gather(*[_sites(p) for p in pdb_ids])
    all_sites = [s for group in per_structure for s in group]
    clusters = cluster_pockets(all_sites)
    pockets = [ConsolidatedPocket(**c) for c in clusters]
    if pockets:
        top = pockets[0]
        interp = (f"{len(pockets)} distinct pocket(s) across {len(pdb_ids)} structures. "
                  f"Most recurrent: a {top.representative_site_type} pocket in {top.structure_count} structure(s) "
                  f"(residues {top.residue_union[0]}–{top.residue_union[-1]}) — the dominant, best-validated site.")
    else:
        interp = f"No non-artifact pockets found across {len(pdb_ids)} structures."
    return ConsolidationResult(uniprot_id=uniprot_id, structures_analyzed=len(pdb_ids),
                               pockets=pockets, interpretation=interp).model_dump()
```
Add `ConsolidatedPocket, ConsolidationResult` to imports.

- [ ] **Step 5: Run to verify pass + suite**

Run: `.venv/bin/python -m pytest tests/test_clients.py::test_cluster_pockets -v && .venv/bin/python -m pytest -q`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add src/pocketscout_mcp/models.py src/pocketscout_mcp/server.py tests/
git commit -m "feat: add ConsolidateBindingSites cross-structure pocket map"
```

---

## Task 18: target_briefing prompt + final pass

**Files:**
- Modify: `src/pocketscout_mcp/server.py`, `README.md`

- [ ] **Step 1: Add the prompt**

In `server.py`, alongside `binding_site_assessment`:
```python
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
```

- [ ] **Step 2: Update README tools/prompt section**

Add a row/sentence documenting the two new tools (`CheckKnownVariants`, `ConsolidateBindingSites`) in the Tools table, and document both prompts (`target_briefing` for quick triage, `binding_site_assessment` for the deep workup) in the Orchestration Prompt section.

- [ ] **Step 3: Verify import + full suite**

Run: `.venv/bin/python -c "import pocketscout_mcp.server" && .venv/bin/python -m pytest -q`
Expected: import OK; all tests pass

- [ ] **Step 4: Smoke-test tool registration**

Run:
```bash
.venv/bin/python -c "
import asyncio, pocketscout_mcp.server as s
async def main():
    tools = await s.mcp.get_tools()
    print(sorted(tools))
asyncio.run(main())
"
```
Expected: lists 8 tools including `CheckKnownVariants` and `ConsolidateBindingSites`. (If `get_tools()` differs in the installed FastMCP version, adjust to the available introspection call.)

- [ ] **Step 5: Commit**

```bash
git add src/pocketscout_mcp/server.py README.md
git commit -m "feat: add target_briefing prompt; document new tools"
```

---

## Final verification

- [ ] Run the complete suite: `.venv/bin/python -m pytest -q` — all green.
- [ ] Confirm 8 tools + 2 prompts register (Task 18 Step 4).
- [ ] `git log --oneline` shows the phased commit history on `triage-updates`.

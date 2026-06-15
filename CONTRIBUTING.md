# Contributing to PocketScout MCP

Thanks for your interest in contributing. This guide covers everything you need to get started.

For the full design philosophy and architecture, read [CLAUDE.md](CLAUDE.md) first — it explains *why* the server is structured the way it is, and will save you time.

---

## Dev Setup

**Requires Python 3.11+.** Your system Python may be older (macOS ships with 3.9). Use [uv](https://docs.astral.sh/uv/) to manage the environment:

```bash
git clone https://github.com/Proprius-Labs/pocketscout-mcp.git
cd pocketscout-mcp

# Create and activate the venv (uv handles the Python version)
uv venv --python 3.11
source .venv/bin/activate   # Windows: .venv\Scripts\activate

# Install the package + dev dependencies
uv pip install -e ".[dev]"
```

The virtual environment lives in `.venv/`. If you prefer to invoke pytest directly without activating, use the full path: `.venv/bin/python -m pytest`.

---

## Running Tests

```bash
.venv/bin/python -m pytest
```

Tests use [respx](https://lundberg.github.io/respx/) to mock all HTTP calls against real-fixture responses — no network access required and no API keys needed. The test suite runs against real EGFR (P00533 / 1M17) and KRAS G12C (P01116 / 6OIM) fixture data.

Run with verbose output to see which tests are executing:

```bash
.venv/bin/python -m pytest -v
```

---

## Test Philosophy

Tests are **real-fixture-driven**, not synthetic. Fixtures were captured from live API responses, so they represent the quirks and inconsistencies of the actual APIs (see the API gotchas in CLAUDE.md).

Every new tool or client must include tests for both:

- **Happy path** — the normal response with complete data.
- **Edge cases** — missing data (no ChEMBL entry, no mouse ortholog, apo structure), network timeouts, and artifact filtering. Missing data is scientifically informative, not just an error, so these paths matter as much as the happy path.

When adding fixture data, capture real API responses and store them as respx mocks in `tests/test_clients.py`. Avoid constructing synthetic responses — they tend to miss the inconsistencies that cause real bugs (e.g., `pdbx_entity_nonpoly` being a dict, not a list).

---

## Code Style and Architecture

PocketScout has some conventions that are not typical Python projects.

### Tool descriptions are prompt engineering

Claude reads the tool docstrings and Pydantic `Field` descriptions to decide *when* and *how* to call each tool. Descriptions must encode scientific reasoning, not just label the data. Compare:

```python
# Bad — describes what, not why
residue_positions: list[int] = Field(description="List of residue positions")

# Good — encodes scientific reasoning
residue_positions: list[int] = Field(
    description="Binding site residue positions from GetBindingSites. "
                "Non-conserved positions here mean your mouse efficacy model "
                "may not predict human response."
)
```

If you add a tool or a model field, the description must be accurate enough that Claude uses it correctly without being told to. Treat them as part of the public API.

### Each tool returns raw data + a pre-computed interpretation

Every tool response includes an `interpretation` field containing a plain-English synthesis of the data with scientific context. This is deliberate: it reduces the amount of reasoning Claude has to do from raw data, and it encodes domain expertise that might not be obvious from the numbers alone. Keep interpretations concise but specific — cite residue counts, potency values, or coverage fractions rather than vague qualitative statements.

### Follow existing client/model/tool patterns

- **Clients** (`src/pocketscout_mcp/clients/`) are async httpx wrappers. Each extends `BaseClient` from `base.py` and raises `APIError` on failure with a scientifically useful message.
- **Models** (`src/pocketscout_mcp/models.py`) are Pydantic v2 models. Field descriptions are the MCP schema — write them accordingly.
- **Tools** (`src/pocketscout_mcp/server.py`) are decorated with `@mcp.tool(name="PascalCase", title="Spaced Title")`. Names are PascalCase; do not change existing names.

### Artifact filtering

`ARTIFACT_LIGANDS` in `pdb.py` filters crystallization artifacts (glycerol, PEG, sulfate, etc.) from binding site analysis. If you find a new artifact contaminant that appears in real structures, add it to the list with a comment explaining what it is.

---

## Pull Request Expectations

Before opening a PR:

1. **All tests pass**: ``.venv/bin/python -m pytest``
2. **Tool/field descriptions are scientifically accurate** — not just syntactically correct.
3. **New tools/clients have both happy-path and edge-case tests.**
4. **Existing client/model/tool patterns are followed** (see above).
5. **Error messages are informative** — state what the absence of data means, not just the HTTP code.

Keep PRs focused. If you're adding a new tool, include the client, the model fields, the server wiring, and the tests in one PR so the reviewer can evaluate the full chain.

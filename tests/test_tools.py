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
    result = await server.get_related_structures(uniprot_id="P00533", limit=5)
    assert result["structures"], "expected at least one structure"
    assert all(s["organism"] == "Homo sapiens" for s in result["structures"])

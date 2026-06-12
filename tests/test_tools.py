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
    result = await server.check_conservation(uniprot_id="P00533", residue_positions=[718, 745])
    species = {s["species"] for s in result["species_results"]}
    assert species == {"mouse", "rat", "cynomolgus"}
    assert all(s["conservation_fraction"] == 1.0 for s in result["species_results"])  # identical seqs


@respx.mock
@pytest.mark.asyncio
async def test_check_known_variants_tool():
    from pocketscout_mcp import server
    respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json=load("uniprot_P00533.json"))
    )
    result = await server.check_known_variants(uniprot_id="P00533", residue_positions=[790])
    assert result["positions_checked"] == 1
    assert isinstance(result["variants"], list)


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
    result = await server.characterize_target(uniprot_id="P01116")
    bands = {r["assessment"] for r in result["alphafold_confidence"]}
    assert "high" in bands and "low" in bands  # multiple regions, not one global band

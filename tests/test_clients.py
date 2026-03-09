"""Tests for PocketScout API clients using real response fixtures."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import respx
import httpx

from pocketscout_mcp.clients.uniprot import UniProtClient, parse_target_profile
from pocketscout_mcp.clients.pdb import PDBClient, parse_structure_metadata, ARTIFACT_LIGANDS
from pocketscout_mcp.clients.alphafold import AlphaFoldClient, analyze_confidence
from pocketscout_mcp.clients.chembl import ChEMBLClient, analyze_bioactivities
from pocketscout_mcp.clients.pubmed import PubMedClient
from pocketscout_mcp.clients.base import APIError

FIXTURES = Path(__file__).parent / "fixtures"


def load_fixture(name: str) -> dict | list:
    with open(FIXTURES / name) as f:
        return json.load(f)


# ---- UniProt ----


class TestUniProtParsing:
    """Test UniProt response parsing with real EGFR data."""

    def setup_method(self):
        self.entry = load_fixture("uniprot_P00533.json")

    def test_parse_target_profile_basic_fields(self):
        profile = parse_target_profile(self.entry)
        assert profile["uniprot_id"] == "P00533"
        assert profile["gene_name"] == "EGFR"
        assert profile["organism"] == "Homo sapiens"
        assert profile["sequence_length"] == 1210

    def test_parse_target_profile_classification(self):
        profile = parse_target_profile(self.entry)
        assert profile["target_class"] == "kinase"

    def test_parse_target_profile_function(self):
        profile = parse_target_profile(self.entry)
        assert "tyrosine kinase" in profile["function_summary"].lower()

    def test_parse_target_profile_location(self):
        profile = parse_target_profile(self.entry)
        assert "membrane" in profile["subcellular_location"].lower()

    def test_parse_target_profile_diseases(self):
        profile = parse_target_profile(self.entry)
        assert len(profile["disease_associations"]) > 0
        assert "Lung cancer" in profile["disease_associations"]

    def test_parse_target_profile_protein_family(self):
        profile = parse_target_profile(self.entry)
        assert "kinase" in profile["protein_family"].lower()


class TestUniProtKRAS:
    """Test UniProt parsing with KRAS data."""

    def setup_method(self):
        self.entry = load_fixture("uniprot_P01116.json")

    def test_kras_basic(self):
        profile = parse_target_profile(self.entry)
        assert profile["gene_name"] == "KRAS"
        assert profile["organism"] == "Homo sapiens"


@respx.mock
@pytest.mark.asyncio
async def test_uniprot_get_entry():
    """Test UniProt client with mocked response."""
    fixture = load_fixture("uniprot_P00533.json")
    respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json=fixture)
    )
    client = UniProtClient()
    entry = await client.get_entry("P00533")
    assert entry["primaryAccession"] == "P00533"
    await client.close()


@respx.mock
@pytest.mark.asyncio
async def test_uniprot_get_sequence():
    """Test sequence extraction."""
    fixture = load_fixture("uniprot_P00533.json")
    respx.get("https://rest.uniprot.org/uniprotkb/P00533.json").mock(
        return_value=httpx.Response(200, json=fixture)
    )
    client = UniProtClient()
    seq = await client.get_sequence("P00533")
    assert len(seq) == 1210
    assert seq[0] == "M"  # EGFR starts with Met
    await client.close()


@respx.mock
@pytest.mark.asyncio
async def test_uniprot_404():
    """Test handling of non-existent entry."""
    respx.get("https://rest.uniprot.org/uniprotkb/FAKE123.json").mock(
        return_value=httpx.Response(404)
    )
    client = UniProtClient()
    with pytest.raises(APIError):
        await client.get_entry("FAKE123")
    await client.close()


# ---- PDB ----


class TestPDBParsing:
    """Test PDB response parsing with real 1M17 data."""

    def setup_method(self):
        self.entry = load_fixture("pdb_1M17_entry.json")

    def test_parse_metadata(self):
        meta = parse_structure_metadata(self.entry)
        assert meta["resolution"] == 2.6
        assert meta["method"] == "X-RAY DIFFRACTION"
        assert "erlotinib" in meta["title"].lower() or "quinazoline" in meta["title"].lower()

    def test_parse_release_date(self):
        meta = parse_structure_metadata(self.entry)
        assert meta["release_date"].startswith("200")  # Released in 2002


class TestPDBNonpolymer:
    """Test nonpolymer entity parsing."""

    def setup_method(self):
        self.entities = load_fixture("pdb_1M17_nonpolymer.json")

    def test_nonpolymer_has_aq4(self):
        comp_ids = []
        for e in self.entities:
            ids = e.get("rcsb_nonpolymer_entity_container_identifiers", {})
            comp_id = ids.get("nonpolymer_comp_id", "")
            if comp_id:
                comp_ids.append(comp_id)
        assert "AQ4" in comp_ids


class TestArtifactFilter:
    """Test the artifact ligand filter list."""

    def test_common_artifacts_filtered(self):
        assert "GOL" in ARTIFACT_LIGANDS  # glycerol
        assert "SO4" in ARTIFACT_LIGANDS  # sulfate
        assert "PEG" in ARTIFACT_LIGANDS  # PEG
        assert "HOH" in ARTIFACT_LIGANDS  # water

    def test_real_drugs_not_filtered(self):
        assert "AQ4" not in ARTIFACT_LIGANDS  # erlotinib
        assert "STI" not in ARTIFACT_LIGANDS  # imatinib
        assert "MOV" not in ARTIFACT_LIGANDS  # sotorasib analog


# ---- AlphaFold ----


class TestAlphaFoldParsing:
    """Test AlphaFold response parsing."""

    def setup_method(self):
        self.pred = load_fixture("alphafold_P00533.json")

    def test_analyze_confidence(self):
        result = analyze_confidence(self.pred, 1210)
        assert result["overall_confidence"] is not None
        assert 70 < result["overall_confidence"] < 90  # EGFR: moderate
        assert len(result["regions"]) > 0
        assert result["regions"][0]["assessment"] == "moderate"

    def test_high_confidence_target(self):
        kras_pred = load_fixture("alphafold_P01116.json")
        result = analyze_confidence(kras_pred, 189)
        assert result["overall_confidence"] > 90  # KRAS: high confidence


# ---- ChEMBL ----


class TestChEMBLAnalysis:
    """Test ChEMBL bioactivity analysis."""

    def setup_method(self):
        self.activities = load_fixture("chembl_P00533_activities.json")
        self.mechanisms = load_fixture("chembl_P00533_mechanisms.json")

    def test_analyze_bioactivities(self):
        result = analyze_bioactivities(self.activities, self.mechanisms)
        assert result["total_compounds_tested"] > 0
        assert result["total_assays"] > 0
        assert result["competitive_landscape"] in (
            "crowded", "moderate", "emerging", "untargeted"
        )

    def test_potency_normalization(self):
        # Activities with nM values should be normalized
        result = analyze_bioactivities(self.activities, [])
        if result["best_potency_nm"] is not None:
            assert result["best_potency_nm"] > 0

    def test_interpretation_generated(self):
        result = analyze_bioactivities(self.activities, self.mechanisms)
        assert len(result["interpretation"]) > 0
        assert "compounds" in result["interpretation"]


class TestChEMBLEmpty:
    """Test ChEMBL analysis with empty data."""

    def test_empty_activities(self):
        result = analyze_bioactivities([], [])
        assert result["competitive_landscape"] == "untargeted"
        assert result["total_compounds_tested"] == 0
        assert result["best_potency_nm"] is None


# ---- PubMed ----


class TestPubMedParsing:
    """Test PubMed response parsing."""

    def setup_method(self):
        self.result = load_fixture("pubmed_EGFR.json")

    def test_papers_parsed(self):
        papers = self.result.get("papers", [])
        assert len(papers) > 0

    def test_paper_fields(self):
        paper = self.result["papers"][0]
        assert "pmid" in paper
        assert "title" in paper
        assert len(paper["title"]) > 0
        assert "authors" in paper

    def test_total_count(self):
        assert self.result["total_found"] > 100  # EGFR: well-studied


# ---- Integration: Conservation helper ----


class TestConservationHelper:
    """Test the context-matching residue finder."""

    def test_find_ortholog_residue_identical_sequences(self):
        from pocketscout_mcp.server import _find_ortholog_residue
        seq = "MADEKVLR"
        result = _find_ortholog_residue(seq, seq, 3)
        assert result == "E"

    def test_find_ortholog_residue_with_insertion(self):
        from pocketscout_mcp.server import _find_ortholog_residue
        human = "MADEKVLRST"
        mouse = "MADXXEKVLRST"  # 2-residue insertion
        # Position 3 in human (E) should map to position 5 in mouse (E)
        result = _find_ortholog_residue(human, mouse, 3)
        assert result == "E"


# ---- Integration: Site classification ----


class TestSiteClassification:
    """Test binding site classification logic."""

    def test_cofactor_classification(self):
        from pocketscout_mcp.server import _classify_site_type
        assert _classify_site_type("ATP", "") == "cofactor"
        assert _classify_site_type("NAD", "") == "cofactor"
        assert _classify_site_type("GDP", "") == "cofactor"

    def test_orthosteric_default(self):
        from pocketscout_mcp.server import _classify_site_type
        assert _classify_site_type("AQ4", "erlotinib") == "orthosteric"

    def test_allosteric_by_name(self):
        from pocketscout_mcp.server import _classify_site_type
        assert _classify_site_type("X01", "allosteric inhibitor") == "allosteric"


class TestConservativeSubstitution:
    """Test conservative substitution checking."""

    def test_identical(self):
        from pocketscout_mcp.server import _is_conservative_substitution
        assert _is_conservative_substitution("L", "L") is True

    def test_conservative(self):
        from pocketscout_mcp.server import _is_conservative_substitution
        assert _is_conservative_substitution("L", "I") is True  # both hydrophobic
        assert _is_conservative_substitution("D", "E") is True  # both negative

    def test_non_conservative(self):
        from pocketscout_mcp.server import _is_conservative_substitution
        assert _is_conservative_substitution("D", "K") is False  # neg vs pos
        assert _is_conservative_substitution("G", "W") is False  # small vs aromatic

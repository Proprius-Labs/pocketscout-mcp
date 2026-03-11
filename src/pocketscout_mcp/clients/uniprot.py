"""UniProt REST API client for protein annotation and sequence data."""

from __future__ import annotations

from .base import BaseClient, APIError


class UniProtClient(BaseClient):
    """Client for the UniProt REST API (2022+ JSON format)."""

    def __init__(self):
        super().__init__(base_url="https://rest.uniprot.org")

    async def get_entry(self, accession: str) -> dict:
        """Fetch a full UniProt entry as JSON."""
        resp = await self.get(f"/uniprotkb/{accession}.json")
        return resp.json()

    async def get_sequence(self, accession: str) -> str:
        """Fetch the amino acid sequence for a UniProt accession."""
        entry = await self.get_entry(accession)
        seq = entry.get("sequence", {})
        if isinstance(seq, dict):
            return seq.get("value", "")
        return ""

    async def get_mouse_ortholog(self, human_accession: str) -> str | None:
        """Find the mouse ortholog for a human UniProt accession.

        Uses UniProt's search to find the equivalent mouse protein.
        Returns the mouse UniProt accession or None.
        """
        # Get gene name from the human entry
        entry = await self.get_entry(human_accession)
        genes = entry.get("genes", [])
        if not genes:
            return None

        gene_name = None
        for gene in genes:
            primary = gene.get("geneName", {})
            if primary:
                gene_name = primary.get("value")
                break

        if not gene_name:
            return None

        # Search for mouse ortholog by gene name
        resp = await self.get(
            "/uniprotkb/search",
            params={
                "query": f"gene_exact:{gene_name} AND organism_id:10090 AND reviewed:true",
                "format": "json",
                "size": "1",
                "fields": "accession",
            },
        )
        data = resp.json()
        results = data.get("results", [])
        if results:
            return results[0].get("primaryAccession")
        return None


def parse_target_profile(entry: dict) -> dict:
    """Extract a structured target profile from a UniProt JSON entry.

    Returns a dict ready to be unpacked into TargetProfile (minus
    AlphaFold fields, which are added separately).
    """
    # Basic identifiers
    accession = entry.get("primaryAccession", "")

    # Gene name
    gene_name = ""
    genes = entry.get("genes", [])
    if genes:
        primary = genes[0].get("geneName", {})
        gene_name = primary.get("value", "") if primary else ""

    # Protein name
    protein_name = ""
    protein_desc = entry.get("proteinDescription", {})
    rec_name = protein_desc.get("recommendedName", {})
    if rec_name:
        full_name = rec_name.get("fullName", {})
        protein_name = full_name.get("value", "") if isinstance(full_name, dict) else ""
    if not protein_name:
        sub_names = protein_desc.get("submittedName", [])
        if sub_names:
            full_name = sub_names[0].get("fullName", {})
            protein_name = full_name.get("value", "") if isinstance(full_name, dict) else ""

    # Organism
    organism = entry.get("organism", {}).get("scientificName", "unknown")

    # Function summary from comments
    function_summary = ""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                function_summary = texts[0].get("value", "")
                break

    # Protein family from comments
    protein_family = ""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SIMILARITY":
            texts = comment.get("texts", [])
            if texts:
                protein_family = texts[0].get("value", "")
                break

    # Subcellular location
    subcellular_location = ""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            locations = comment.get("subcellularLocations", [])
            if locations:
                loc_parts = []
                for loc in locations:
                    loc_val = loc.get("location", {})
                    if isinstance(loc_val, dict):
                        val = loc_val.get("value", "")
                        if val:
                            loc_parts.append(val)
                subcellular_location = "; ".join(loc_parts)
            break

    # Disease associations
    disease_associations = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "DISEASE":
            disease = comment.get("disease", {})
            if disease:
                disease_name = disease.get("diseaseId", "")
                if disease_name:
                    disease_associations.append(disease_name)

    # Sequence length
    seq_info = entry.get("sequence", {})
    sequence_length = seq_info.get("length", 0) if isinstance(seq_info, dict) else 0

    # Topology — extracellular, transmembrane, cytoplasmic regions
    topology = _parse_topology(entry)

    # Classify target
    target_class = _classify_target(protein_family, protein_name, entry)

    return {
        "uniprot_id": accession,
        "gene_name": gene_name,
        "protein_name": protein_name,
        "organism": organism,
        "function_summary": function_summary,
        "protein_family": protein_family,
        "subcellular_location": subcellular_location,
        "disease_associations": disease_associations,
        "target_class": target_class,
        "sequence_length": sequence_length,
        "topology": topology,
    }


def _parse_topology(entry: dict) -> list[dict]:
    """Extract topological regions from UniProt features.

    Identifies extracellular, transmembrane, and cytoplasmic regions.
    Returns a list of dicts ready to be unpacked into TopologyRegion.
    """
    type_map = {
        "Topological domain": True,
        "Transmembrane": True,
        "Signal": True,
    }

    desc_to_region = {
        "extracellular": "extracellular",
        "cytoplasmic": "cytoplasmic",
        "lumenal": "extracellular",  # ER lumen is topologically extracellular
        "periplasmic": "extracellular",
        "helical": "transmembrane",
    }

    regions = []
    for feature in entry.get("features", []):
        feat_type = feature.get("type", "")
        if feat_type not in type_map:
            continue

        loc = feature.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue

        desc = feature.get("description", "").lower().strip()

        if feat_type == "Signal":
            region_type = "signal_peptide"
        elif feat_type == "Transmembrane":
            region_type = "transmembrane"
        else:
            region_type = desc_to_region.get(desc, desc or "unknown")

        regions.append({
            "start": start,
            "end": end,
            "region_type": region_type,
            "description": feature.get("description", ""),
        })

    return regions


def _classify_target(family: str, name: str, entry: dict) -> str:
    """Classify the target into a druggability class."""
    text = f"{family} {name}".lower()

    if "kinase" in text:
        return "kinase"
    if "g protein-coupled receptor" in text or "gpcr" in text:
        return "gpcr"
    if "protease" in text or "peptidase" in text:
        return "protease"
    if "ion channel" in text or "voltage-gated" in text:
        return "ion_channel"
    if "nuclear receptor" in text:
        return "nuclear_receptor"

    # Check for enzyme-related keywords
    if any(kw in text for kw in ("enzyme", "transferase", "hydrolase", "oxidoreductase", "ligase", "lyase", "isomerase", "synthase", "synthetase")):
        return "enzyme_other"

    # Check keywords for PPI targets
    keywords = entry.get("keywords", [])
    kw_values = [kw.get("name", "").lower() for kw in keywords]
    if any("receptor" in kw for kw in kw_values):
        if any("cell surface" in kw or "membrane" in kw for kw in kw_values):
            return "ppi"

    return "other"

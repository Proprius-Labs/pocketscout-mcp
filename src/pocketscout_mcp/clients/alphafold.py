"""AlphaFold Database API client for structure confidence data."""

from __future__ import annotations

from .base import BaseClient, APIError


class AlphaFoldClient(BaseClient):
    """Client for the AlphaFold DB REST API."""

    def __init__(self):
        super().__init__(
            base_url="https://alphafold.ebi.ac.uk",
            timeout=15.0,
        )

    async def get_prediction(self, uniprot_id: str) -> dict:
        """Fetch AlphaFold prediction metadata for a UniProt accession.

        The API returns a list — we take the first element (latest model).
        """
        resp = await self.get(f"/api/prediction/{uniprot_id}")
        data = resp.json()
        # API returns a list; take the first (latest) prediction
        if isinstance(data, list):
            if not data:
                raise APIError(f"No AlphaFold prediction for {uniprot_id}")
            return data[0]
        return data


def analyze_confidence(prediction: dict, sequence_length: int) -> dict:
    """Analyze AlphaFold confidence from prediction metadata.

    Returns overall confidence and region-level assessments.
    The global pLDDT is available in the prediction metadata;
    per-residue pLDDT would require downloading and parsing the CIF file.
    """
    # Global pLDDT is in the prediction metadata
    global_plddt = prediction.get("globalMetricValue")
    if global_plddt is None:
        # Try alternative field names
        global_plddt = prediction.get("plddt")

    # Without per-residue data from the CIF, we create a single-region
    # assessment based on the global metric
    regions = []
    warnings = []

    if global_plddt is not None:
        if global_plddt >= 90:
            assessment = "high"
        elif global_plddt >= 70:
            assessment = "moderate"
        else:
            assessment = "low"

        regions.append({
            "start": 1,
            "end": sequence_length,
            "mean_plddt": round(global_plddt, 1),
            "assessment": assessment,
        })

        if global_plddt < 70:
            warnings.append(
                f"Low overall AlphaFold confidence ({global_plddt:.0f} pLDDT). "
                "Structure predictions may be unreliable — prioritize experimental structures."
            )
        elif global_plddt < 80:
            warnings.append(
                f"Moderate AlphaFold confidence ({global_plddt:.0f} pLDDT). "
                "Some regions may have uncertain structure. Cross-reference with experimental data."
            )
    else:
        warnings.append("Global pLDDT score not available in AlphaFold metadata.")

    return {
        "overall_confidence": global_plddt,
        "regions": regions,
        "warnings": warnings,
    }

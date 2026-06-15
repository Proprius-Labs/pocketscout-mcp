"""AlphaFold Database API client for structure confidence data."""

from __future__ import annotations

from .base import BaseClient, APIError


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
        data = await self.get_json(f"/api/prediction/{uniprot_id}")
        # API returns a list; take the first (latest) prediction
        if isinstance(data, list):
            if not data:
                raise APIError(f"No AlphaFold prediction for {uniprot_id}")
            return data[0]
        return data

    async def get_per_residue_plddt(self, uniprot_id: str) -> list[float]:
        """Download the AlphaFold CIF and read per-residue pLDDT (B-factor column).

        Best-effort: returns [] on any failure (no prediction, no cifUrl,
        download error, or parse error) so callers can fall back to the
        global confidence path.
        """
        import httpx as _httpx
        try:
            prediction = await self.get_prediction(uniprot_id)
            cif_url = prediction.get("cifUrl")
            if not cif_url:
                return []
            async with _httpx.AsyncClient(timeout=30.0, follow_redirects=True) as dl:
                resp = await dl.get(cif_url)
            if resp.status_code != 200:
                return []
            import gemmi
            st = gemmi.make_structure_from_block(gemmi.cif.read_string(resp.text).sole_block())
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


def analyze_confidence(prediction: dict, sequence_length: int, per_residue: list[float] | None = None) -> dict:
    """Analyze AlphaFold confidence from prediction metadata.

    Returns overall confidence and region-level assessments.
    When per_residue values are supplied (from the CIF B-factor column), produces
    per-band segmentation with region-level low-confidence warnings.
    Falls back to a single-region assessment using the global pLDDT metric.
    """
    # Per-residue path: segment by confidence band, build region-level warnings
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

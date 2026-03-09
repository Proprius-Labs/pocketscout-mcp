"""ChEMBL API client for bioactivity and competitive landscape data."""

from __future__ import annotations

import statistics

from .base import BaseClient, APIError


class ChEMBLClient(BaseClient):
    """Client for the ChEMBL REST API.

    ChEMBL can be slow for heavily studied targets (30-60s).
    Uses longer timeouts.
    """

    def __init__(self):
        super().__init__(
            base_url="https://www.ebi.ac.uk/chembl/api/data",
            timeout=60.0,
        )

    async def get_target_by_uniprot(self, uniprot_id: str) -> dict | None:
        """Map a UniProt accession to a ChEMBL target.

        Returns the first matching ChEMBL target record, or None.
        """
        try:
            resp = await self.get(
                "/target.json",
                params={
                    "target_components__accession": uniprot_id,
                    "limit": "1",
                },
            )
        except APIError:
            return None

        data = resp.json()
        targets = data.get("targets", [])
        return targets[0] if targets else None

    async def get_bioactivities(
        self,
        chembl_target_id: str,
        limit: int = 500,
    ) -> list[dict]:
        """Fetch bioactivity data (IC50, Ki, Kd, EC50) for a target.

        Returns up to `limit` activity records. For heavily studied
        targets, this may be a subset of all available data.
        """
        try:
            resp = await self.get(
                "/activity.json",
                params={
                    "target_chembl_id": chembl_target_id,
                    "standard_type__in": "IC50,Ki,Kd,EC50",
                    "limit": str(limit),
                },
            )
        except APIError:
            return []

        data = resp.json()
        return data.get("activities", [])

    async def get_clinical_candidates(self, chembl_target_id: str) -> list[dict]:
        """Fetch mechanism-of-action records (clinical/approved drugs).

        Enriches with molecule preferred names via batch lookup.
        """
        try:
            resp = await self.get(
                "/mechanism.json",
                params={
                    "target_chembl_id": chembl_target_id,
                    "limit": "100",
                },
            )
        except APIError:
            return []

        data = resp.json()
        mechanisms = data.get("mechanisms", [])

        # Collect unique molecule IDs for name resolution
        unique_mol_ids = list({
            m.get("molecule_chembl_id", "")
            for m in mechanisms
            if m.get("molecule_chembl_id")
        })

        # Batch resolve names (limit to first 20 to avoid excessive API calls)
        name_map: dict[str, str] = {}
        for mol_id in unique_mol_ids[:20]:
            try:
                mol_resp = await self.get(f"/molecule/{mol_id}.json")
                mol_data = mol_resp.json()
                pref_name = mol_data.get("pref_name")
                if pref_name:
                    name_map[mol_id] = pref_name
            except APIError:
                pass

        # Attach resolved names
        for mech in mechanisms:
            mol_id = mech.get("molecule_chembl_id", "")
            if mol_id in name_map:
                mech["_resolved_name"] = name_map[mol_id]

        return mechanisms


def analyze_bioactivities(activities: list[dict], mechanisms: list[dict]) -> dict:
    """Analyze bioactivity data to characterize competitive landscape.

    Returns a dict ready to be unpacked into LigandHistory (minus
    uniprot_id and chembl_target_id).
    """
    # Collect potency values, normalizing to nM
    potencies_nm: list[float] = []
    compound_ids: set[str] = set()
    assay_ids: set[str] = set()

    for act in activities:
        mol_id = act.get("molecule_chembl_id", "")
        if mol_id:
            compound_ids.add(mol_id)

        assay_id = act.get("assay_chembl_id", "")
        if assay_id:
            assay_ids.add(assay_id)

        # Get potency value
        value = act.get("standard_value")
        units = act.get("standard_units", "")

        if value is not None:
            try:
                val = float(value)
            except (ValueError, TypeError):
                continue

            # Normalize to nM
            if units == "nM":
                potencies_nm.append(val)
            elif units == "uM":
                potencies_nm.append(val * 1000)
            elif units == "pM":
                potencies_nm.append(val / 1000)
            elif units == "M":
                potencies_nm.append(val * 1e9)
            elif units == "mM":
                potencies_nm.append(val * 1e6)

    # Clinical candidates — prefer resolved drug names over ChEMBL IDs
    clinical_names = []
    seen_mol_ids: set[str] = set()
    for mech in mechanisms:
        mol_id = mech.get("molecule_chembl_id", "")
        if mol_id in seen_mol_ids:
            continue
        seen_mol_ids.add(mol_id)
        # Prefer _resolved_name (from batch lookup), then molecule_pref_name, then ChEMBL ID
        pref_name = (
            mech.get("_resolved_name")
            or mech.get("molecule_pref_name")
        )
        if pref_name and pref_name not in clinical_names:
            clinical_names.append(pref_name)

    # Calculate metrics
    total_compounds = len(compound_ids)
    total_assays = len(assay_ids)
    best_potency = min(potencies_nm) if potencies_nm else None
    median_potency = statistics.median(potencies_nm) if potencies_nm else None

    # Classify landscape
    # Count compounds with submicromolar activity
    potent_count = sum(1 for p in potencies_nm if p < 1000)
    if potent_count > 100:
        landscape = "crowded"
    elif potent_count > 10:
        landscape = "moderate"
    elif potent_count > 0:
        landscape = "emerging"
    else:
        landscape = "untargeted"

    # Interpretation
    interp_parts = [f"{total_compounds} compounds tested across {total_assays} assays."]

    if best_potency is not None:
        if best_potency < 10:
            interp_parts.append(f"Best potency: {best_potency:.1f} nM (very potent — sub-10 nM leads exist).")
        elif best_potency < 100:
            interp_parts.append(f"Best potency: {best_potency:.1f} nM (potent leads available).")
        elif best_potency < 1000:
            interp_parts.append(f"Best potency: {best_potency:.0f} nM (moderate potency — room for optimization).")
        else:
            interp_parts.append(f"Best potency: {best_potency:.0f} nM (weak — significant optimization needed).")

    if clinical_names:
        interp_parts.append(f"Clinical candidates: {', '.join(clinical_names[:5])}.")

    landscape_desc = {
        "crowded": "Crowded landscape — many potent compounds exist. De novo design should target novel binding modes, allosteric sites, or differentiated modalities.",
        "moderate": "Moderate competitive landscape. Known pharmacology exists but opportunity remains for improved selectivity or novel mechanisms.",
        "emerging": "Emerging landscape — limited but active pharmacology. Good opportunity for first-in-class approaches.",
        "untargeted": "No significant bioactivity data. Novel target — high opportunity but less prior validation.",
    }
    interp_parts.append(landscape_desc[landscape])

    return {
        "total_compounds_tested": total_compounds,
        "total_assays": total_assays,
        "best_potency_nm": round(best_potency, 2) if best_potency is not None else None,
        "median_potency_nm": round(median_potency, 2) if median_potency is not None else None,
        "clinical_candidates": clinical_names,
        "competitive_landscape": landscape,
        "interpretation": " ".join(interp_parts),
    }

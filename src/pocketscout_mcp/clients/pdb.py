"""RCSB PDB API client for structure and binding site data."""

from __future__ import annotations

from .base import BaseClient, APIError


# Crystallization artifacts — these are NOT real binding events
ARTIFACT_LIGANDS = {
    "GOL", "EDO", "PEG", "PG4", "1PE", "P6G", "MPD",  # Polyethylene glycols, glycerol
    "SO4", "PO4", "CIT", "ACT", "FMT", "OXL",          # Buffer/precipitant ions
    "DMS", "MSE",                                         # DMSO, selenomethionine
    "BME", "EOH", "MES", "TRS",                          # Buffers
    "CL", "NA", "K", "MG", "ZN", "CA", "MN", "FE",     # Common ions
    "IOD", "BR",                                          # Halides from soaking
    "HOH", "WAT", "DOD",                                  # Water
    "UNX", "UNL",                                         # Unknown atoms/ligands
}


def parse_structure_metadata(entry: dict) -> dict:
    """Extract key metadata from a PDB entry response.

    Returns a dict with title, resolution, method, and release_date.
    """
    title = entry.get("struct", {}).get("title", "")

    # Resolution — try multiple paths
    resolution = None
    refine = entry.get("refine", [])
    if isinstance(refine, list) and refine:
        resolution = refine[0].get("ls_d_res_high")
    if resolution is None:
        resolution = entry.get("rcsb_entry_info", {}).get("resolution_combined", [None])
        if isinstance(resolution, list) and resolution:
            resolution = resolution[0]
        elif isinstance(resolution, (int, float)):
            pass
        else:
            resolution = None
    if resolution is not None:
        try:
            resolution = float(resolution)
        except (ValueError, TypeError):
            resolution = None

    # Method
    method = ""
    exptl = entry.get("exptl", [])
    if isinstance(exptl, list) and exptl:
        method = exptl[0].get("method", "")

    # Release date
    release_date = entry.get("rcsb_accession_info", {}).get("initial_release_date", "")
    if release_date and "T" in release_date:
        release_date = release_date.split("T")[0]

    return {
        "title": title,
        "resolution": resolution,
        "method": method,
        "release_date": release_date,
    }


class PDBClient(BaseClient):
    """Client for the RCSB PDB Data and Search APIs."""

    def __init__(self):
        super().__init__(
            base_url="https://data.rcsb.org",
            timeout=30.0,
        )

    async def get_entry(self, pdb_id: str) -> dict:
        """Fetch core entry metadata for a PDB structure."""
        resp = await self.get(f"/rest/v1/core/entry/{pdb_id.upper()}")
        return resp.json()

    async def get_uniprot_residue_range(self, pdb_id: str) -> tuple[int, int] | None:
        """Get the UniProt residue range covered by a PDB structure.

        Returns (start, end) tuple in UniProt numbering, or None if
        the mapping can't be determined. Uses the polymer entity
        alignment data from RCSB.
        """
        pdb_id = pdb_id.upper()
        try:
            entry = await self.get_entry(pdb_id)
        except APIError:
            return None

        entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
            "polymer_entity_ids", []
        )

        for eid in entity_ids:
            try:
                resp = await self.get(f"/rest/v1/core/polymer_entity/{pdb_id}/{eid}")
                entity = resp.json()

                # Check if this entity maps to UniProt
                identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
                if not identifiers.get("uniprot_ids"):
                    continue

                # Get alignment info from reference sequence identifiers
                align = entity.get("rcsb_polymer_entity_align", [])
                for a in align:
                    if a.get("reference_database_name") == "UniProt":
                        aligned = a.get("aligned_regions", [])
                        if aligned:
                            starts = [r.get("ref_beg_seq_id", 0) for r in aligned]
                            ends = [r.get("ref_beg_seq_id", 0) + r.get("length", 0) - 1 for r in aligned]
                            return (min(starts), max(ends))
            except APIError:
                continue

        return None

    async def get_uniprot_mapping(self, pdb_id: str) -> list[str]:
        """Get UniProt accessions mapped to a PDB structure.

        Uses the polymer_entities endpoint to find UniProt references.
        """
        pdb_id = pdb_id.upper()
        # First get entity IDs
        try:
            entry = await self.get_entry(pdb_id)
        except APIError:
            return []

        entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        accessions = []

        for eid in entity_ids:
            try:
                resp = await self.get(f"/rest/v1/core/polymer_entity/{pdb_id}/{eid}")
                entity = resp.json()
                # UniProt references are in rcsb_polymer_entity_container_identifiers
                identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
                uniprot_ids = identifiers.get("uniprot_ids", [])
                accessions.extend(uniprot_ids)
            except APIError:
                continue

        return list(dict.fromkeys(accessions))  # deduplicate preserving order

    async def search_by_uniprot(self, uniprot_id: str, limit: int = 20) -> list[str]:
        """Search PDB for all structures of a UniProt accession.

        Uses the RCSB Search API v2.
        """
        import httpx

        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id,
                },
            },
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": limit,
                },
                "results_content_type": ["experimental"],
                "sort": [
                    {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
                ],
            },
            "return_type": "entry",
        }

        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.post(
                "https://search.rcsb.org/rcsbsearch/v2/query",
                json=query,
            )
            if resp.status_code == 204:
                # No results
                return []
            resp.raise_for_status()
            data = resp.json()

        results = data.get("result_set", [])
        return [r.get("identifier", "") for r in results if r.get("identifier")]

    async def get_nonpolymer_entities(self, pdb_id: str) -> list[dict]:
        """Get all non-polymer entities (ligands) for a structure."""
        pdb_id = pdb_id.upper()
        try:
            entry = await self.get_entry(pdb_id)
        except APIError:
            return []

        entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
            "non_polymer_entity_ids", []
        )
        entities = []
        for eid in entity_ids:
            try:
                resp = await self.get(f"/rest/v1/core/nonpolymer_entity/{pdb_id}/{eid}")
                entities.append(resp.json())
            except APIError:
                continue
        return entities

    async def get_binding_sites(self, pdb_id: str) -> dict:
        """Get binding site information with residue contacts from coordinates.

        Downloads the mmCIF file and uses gemmi to compute which protein
        residues are within contact distance (4.5 A) of each non-artifact
        ligand. Returns ligand info, contact residues, and filtered artifacts.
        """
        pdb_id = pdb_id.upper()
        entities = await self.get_nonpolymer_entities(pdb_id)

        ligands = []
        artifacts_filtered = []

        for entity in entities:
            comp_id = ""
            identifiers = entity.get("rcsb_nonpolymer_entity_container_identifiers", {})
            comp_id = identifiers.get("nonpolymer_comp_id", "") or identifiers.get("comp_id", "")

            if not comp_id:
                pdbx = entity.get("pdbx_entity_nonpoly", {})
                if isinstance(pdbx, list) and pdbx:
                    comp_id = pdbx[0].get("comp_id", "")
                elif isinstance(pdbx, dict):
                    comp_id = pdbx.get("comp_id", "")

            if not comp_id:
                continue

            if comp_id in ARTIFACT_LIGANDS:
                artifacts_filtered.append(comp_id)
                continue

            name = ""
            nonpoly = entity.get("pdbx_entity_nonpoly", {})
            if isinstance(nonpoly, list) and nonpoly:
                name = nonpoly[0].get("name", "")
            elif isinstance(nonpoly, dict):
                name = nonpoly.get("name", "")
            if not name:
                rcsb = entity.get("rcsb_nonpolymer_entity", {})
                name = rcsb.get("pdbx_description", "")

            ligands.append({
                "comp_id": comp_id,
                "name": name,
            })

        # Compute residue contacts from coordinates
        contacts_by_ligand = await self._compute_residue_contacts(pdb_id, ligands)

        # Attach contacts to ligands
        for lig in ligands:
            lig["contacts"] = contacts_by_ligand.get(lig["comp_id"], [])

        return {
            "ligands": ligands,
            "artifacts_filtered": artifacts_filtered,
        }

    async def _compute_residue_contacts(
        self,
        pdb_id: str,
        ligands: list[dict],
        cutoff: float = 4.5,
    ) -> dict[str, list[dict]]:
        """Download mmCIF and compute protein residue contacts for each ligand.

        Returns a dict mapping comp_id to a list of contacting residues:
        [{"chain": "A", "residue_name": "LEU", "residue_number": 694}, ...]
        """
        if not ligands:
            return {}

        try:
            import gemmi
        except ImportError:
            # gemmi not installed — return empty contacts
            return {}

        ligand_ids = {lig["comp_id"] for lig in ligands}

        # Standard amino acid residue names
        amino_acids = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL", "MSE",  # selenomethionine
        }

        try:
            import httpx as _httpx
            async with _httpx.AsyncClient(timeout=30.0) as dl_client:
                resp = await dl_client.get(
                    f"https://files.rcsb.org/download/{pdb_id}.cif"
                )
                if resp.status_code != 200:
                    return {}
                cif_text = resp.text
        except Exception:
            return {}

        try:
            doc = gemmi.cif.read_string(cif_text)
            block = doc.sole_block()
            st = gemmi.make_structure_from_block(block)
            model = st[0]

            ns = gemmi.NeighborSearch(model, st.cell, cutoff)
            ns.populate(include_h=False)

            contacts_by_ligand: dict[str, list[dict]] = {}

            for chain in model:
                for residue in chain:
                    if residue.name not in ligand_ids:
                        continue
                    comp_id = residue.name
                    seen = set()
                    for atom in residue:
                        for mark in ns.find_atoms(atom.pos, "\0", radius=cutoff):
                            cra = mark.to_cra(model)
                            if cra.residue.name not in amino_acids:
                                continue
                            key = (cra.chain.name, cra.residue.name, str(cra.residue.seqid))
                            if key in seen:
                                continue
                            seen.add(key)

                    contact_list = []
                    for ch, resname, seqid in sorted(
                        seen, key=lambda x: (x[0], int(x[2]) if x[2].isdigit() else 0)
                    ):
                        seq_num = int(seqid) if seqid.isdigit() else 0
                        contact_list.append({
                            "chain": ch,
                            "residue_name": resname,
                            "residue_number": seq_num,
                        })

                    if comp_id not in contacts_by_ligand:
                        contacts_by_ligand[comp_id] = contact_list
                    else:
                        # Merge contacts from multiple instances of same ligand
                        existing_keys = {
                            (c["chain"], c["residue_number"])
                            for c in contacts_by_ligand[comp_id]
                        }
                        for c in contact_list:
                            if (c["chain"], c["residue_number"]) not in existing_keys:
                                contacts_by_ligand[comp_id].append(c)

            return contacts_by_ligand

        except Exception:
            return {}

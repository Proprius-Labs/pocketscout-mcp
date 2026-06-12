"""PubMed E-utilities client for literature search."""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET

from .base import BaseClient, APIError


def _parse_abstracts(xml_text: str) -> dict[str, str]:
    """Parse PubMed efetch XML and return {pmid: full_abstract_text}."""
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


class PubMedClient(BaseClient):
    """Client for NCBI PubMed E-utilities API."""

    def __init__(self):
        super().__init__(
            base_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
            timeout=15.0,
        )
        self._api_key = os.environ.get("NCBI_API_KEY")

    def _common_params(self) -> dict:
        """Common parameters for all E-utility requests."""
        params: dict[str, str] = {"retmode": "json"}
        if self._api_key:
            params["api_key"] = self._api_key
        return params

    async def fetch_abstracts(self, pmids: list[str]) -> dict[str, str]:
        """Fetch full abstracts for a list of PMIDs via efetch (XML).

        Returns a dict of {pmid: abstract_text}. Returns empty dict on failure
        so callers degrade gracefully.
        """
        if not pmids:
            return {}
        params = self._common_params()
        params.update({
            "db": "pubmed",
            "id": ",".join(pmids),
            "rettype": "abstract",
            "retmode": "xml",
        })
        try:
            resp = await self.get("/efetch.fcgi", params=params)
        except APIError:
            return {}
        return _parse_abstracts(resp.text)

    async def search_target_literature(
        self,
        gene_name: str,
        context: str | None = None,
        max_results: int = 10,
    ) -> dict:
        """Search PubMed for structural/design papers about a target.

        Builds a focused query combining gene name with structural
        biology and drug design terms.
        """
        # Build search query focused on structural biology and drug design
        base_terms = f'({gene_name}[Title/Abstract])'
        structure_terms = '(binding site[Title/Abstract] OR crystal structure[Title/Abstract] OR "drug design"[Title/Abstract] OR allosteric[Title/Abstract] OR druggable[Title/Abstract] OR "de novo design"[Title/Abstract])'

        query = f"{base_terms} AND {structure_terms}"

        if context:
            query = f"{query} AND ({context}[Title/Abstract])"

        # Search
        params = self._common_params()
        params.update({
            "db": "pubmed",
            "term": query,
            "retmax": str(max_results),
            "sort": "relevance",
        })

        resp = await self.get("/esearch.fcgi", params=params)
        search_data = resp.json()

        esearch_result = search_data.get("esearchresult", {})
        id_list = esearch_result.get("idlist", [])
        total_count = int(esearch_result.get("count", 0))

        if not id_list:
            return {
                "query": query,
                "total_found": total_count,
                "papers": [],
            }

        # Fetch summaries
        summary_params = self._common_params()
        summary_params.update({
            "db": "pubmed",
            "id": ",".join(id_list),
        })

        resp = await self.get("/esummary.fcgi", params=summary_params)
        summary_data = resp.json()

        result = summary_data.get("result", {})
        uids = result.get("uids", [])

        papers = []
        for uid in uids:
            article = result.get(uid, {})
            if not isinstance(article, dict):
                continue

            # Parse authors
            authors_list = article.get("authors", [])
            if authors_list:
                first_author = authors_list[0].get("name", "Unknown")
                author_str = f"{first_author} et al." if len(authors_list) > 1 else first_author
            else:
                author_str = "Unknown"

            # Parse year from pubdate
            pubdate = article.get("pubdate", "")
            year = None
            if pubdate:
                parts = pubdate.split()
                if parts:
                    try:
                        year = int(parts[0])
                    except ValueError:
                        pass

            papers.append({
                "pmid": uid,
                "title": article.get("title", ""),
                "authors": author_str,
                "year": year,
                "journal": article.get("source", ""),
                "abstract_snippet": "",
            })

        abstracts = await self.fetch_abstracts(id_list)
        for p in papers:
            p["abstract_snippet"] = abstracts.get(p["pmid"], "")[:300]

        return {
            "query": query,
            "total_found": total_count,
            "papers": papers,
        }

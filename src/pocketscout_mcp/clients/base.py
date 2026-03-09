"""Shared async HTTP client with retry logic for biological databases."""

from __future__ import annotations

import httpx


class APIError(Exception):
    """Raised when an API request fails after retries."""

    def __init__(self, message: str, status_code: int | None = None):
        self.status_code = status_code
        super().__init__(message)


class BaseClient:
    """Async HTTP client with retries and timeout handling.

    All PocketScout API clients inherit from this. Biological databases
    can be slow and flaky — retries and generous timeouts are essential.
    """

    def __init__(
        self,
        base_url: str,
        timeout: float = 30.0,
        max_retries: int = 3,
    ):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.max_retries = max_retries
        self._client: httpx.AsyncClient | None = None

    async def _get_client(self) -> httpx.AsyncClient:
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(
                base_url=self.base_url,
                timeout=httpx.Timeout(self.timeout, connect=10.0),
                headers={"Accept": "application/json"},
                follow_redirects=True,
            )
        return self._client

    async def _request(
        self,
        method: str,
        path: str,
        **kwargs,
    ) -> httpx.Response:
        client = await self._get_client()
        last_exc: Exception | None = None

        for attempt in range(self.max_retries):
            try:
                resp = await client.request(method, path, **kwargs)
                if resp.status_code == 404:
                    raise APIError(f"Not found: {path}", status_code=404)
                if resp.status_code == 429:
                    # Rate limited — wait briefly and retry
                    import asyncio
                    await asyncio.sleep(1.0 * (attempt + 1))
                    continue
                resp.raise_for_status()
                return resp
            except APIError:
                raise
            except httpx.HTTPStatusError as e:
                last_exc = e
                if e.response.status_code >= 500:
                    import asyncio
                    await asyncio.sleep(1.0 * (attempt + 1))
                    continue
                raise APIError(
                    f"HTTP {e.response.status_code} from {path}",
                    status_code=e.response.status_code,
                ) from e
            except (httpx.TimeoutException, httpx.ConnectError) as e:
                last_exc = e
                import asyncio
                await asyncio.sleep(1.0 * (attempt + 1))
                continue

        raise APIError(f"Request failed after {self.max_retries} retries: {last_exc}")

    async def get(self, path: str, **kwargs) -> httpx.Response:
        return await self._request("GET", path, **kwargs)

    async def post(self, path: str, **kwargs) -> httpx.Response:
        return await self._request("POST", path, **kwargs)

    async def close(self):
        if self._client and not self._client.is_closed:
            await self._client.aclose()

#!/usr/bin/env python3
"""Fetch BRCA1 data from the Ensembl REST API and save JSON + FASTA outputs.

This script looks up the BRCA1 gene symbol in Homo sapiens, retrieves the full
record for the corresponding Ensembl stable ID, and downloads sequence data as
FASTA files. Results are written to an output directory (default: ./outputs).
"""

from __future__ import annotations

import argparse
import json
import pathlib
import ssl
import sys
import time
from typing import Any, Dict, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

try:
    import certifi
except ModuleNotFoundError:  # pragma: no cover - optional dependency
    certifi = None

ENSEMBL_REST_BASE = "https://rest.ensembl.org"
USER_AGENT = "ensembl-brca1-sequence-fetcher/1.0"


class EnsemblClient:
    """Minimal wrapper around the Ensembl REST API."""

    def __init__(
        self,
        base_url: str = ENSEMBL_REST_BASE,
        max_retries: int = 3,
        backoff_seconds: float = 1.0,
        verify_ssl: bool = True,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds
        self.verify_ssl = verify_ssl

    def _request(self, path: str, *, accept: str, params: Optional[Dict[str, Any]] = None) -> bytes:
        url = f"{self.base_url}{path}"
        if params:
            url = f"{url}?{urlencode(params)}"

        request = Request(url, headers={"Accept": accept, "User-Agent": USER_AGENT})
        context = self._build_ssl_context()

        for attempt in range(1, self.max_retries + 1):
            try:
                with urlopen(request, timeout=30, context=context) as response:  # nosec B310
                    return response.read()
            except HTTPError as error:
                if 500 <= error.code < 600 and attempt < self.max_retries:
                    time.sleep(self.backoff_seconds * attempt)
                    continue
                raise
            except URLError:
                if attempt < self.max_retries:
                    time.sleep(self.backoff_seconds * attempt)
                    continue
                raise
        raise RuntimeError("Exceeded maximum retries")

    def _build_ssl_context(self) -> ssl.SSLContext:
        if not self.verify_ssl:
            return ssl._create_unverified_context()  # noqa: SLF001 - intentional opt-out
        if certifi:
            return ssl.create_default_context(cafile=certifi.where())
        return ssl.create_default_context()

    def lookup_symbol(self, species: str, symbol: str) -> Dict[str, Any]:
        payload = self._request(
            f"/lookup/symbol/{species}/{symbol}",
            accept="application/json",
            params={"expand": 1},
        )
        return json.loads(payload)

    def lookup_id(self, stable_id: str) -> Dict[str, Any]:
        payload = self._request(
            f"/lookup/id/{stable_id}",
            accept="application/json",
            params={"expand": 1},
        )
        return json.loads(payload)

    def fetch_sequence_fasta(self, stable_id: str, *, seq_type: str) -> str:
        candidates = [stable_id]
        if "." in stable_id:
            candidates.append(stable_id.split(".")[0])

        last_error: Optional[Exception] = None
        for candidate in candidates:
            try:
                payload = self._request(
                    f"/sequence/id/{candidate}",
                    accept="text/x-fasta",
                    params={"type": seq_type},
                )
                return payload.decode("utf-8")
            except HTTPError as error:
                last_error = error
            except URLError as error:
                last_error = error

        if last_error:
            raise last_error
        raise RuntimeError("Unexpected error fetching sequence")


def write_text(path: pathlib.Path, data: str) -> None:
    path.write_text(data, encoding="utf-8")


def write_json(path: pathlib.Path, data: Dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Fetch BRCA1 data from Ensembl REST API")
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("outputs"),
        help="Directory where fetched files are written (default: ./outputs)",
    )
    parser.add_argument(
        "--symbol",
        default="BRCA1",
        help="Gene symbol to fetch (default: BRCA1)",
    )
    parser.add_argument(
        "--species",
        default="homo_sapiens",
        help="Ensembl species name (default: homo_sapiens)",
    )
    parser.add_argument(
        "--skip-protein",
        action="store_true",
        help="Skip fetching the canonical protein FASTA",
    )
    parser.add_argument(
        "--insecure",
        action="store_true",
        help="Disable TLS certificate verification (not recommended).",
    )
    args = parser.parse_args(argv)

    output_dir = args.output
    output_dir.mkdir(parents=True, exist_ok=True)

    client = EnsemblClient(verify_ssl=not args.insecure)

    try:
        lookup = client.lookup_symbol(args.species, args.symbol)
    except HTTPError as error:
        sys.stderr.write(f"Failed to resolve symbol {args.symbol!r}: HTTP {error.code}\n")
        return 1
    except URLError as error:
        sys.stderr.write(f"Network error while resolving symbol {args.symbol!r}: {error}\n")
        return 1

    stable_id = lookup.get("id")
    if not stable_id:
        sys.stderr.write("No stable ID found in lookup response.\n")
        return 1

    write_json(output_dir / "brca1_lookup.json", lookup)

    try:
        full_record = client.lookup_id(stable_id)
    except HTTPError as error:
        sys.stderr.write(f"Failed to fetch record for {stable_id}: HTTP {error.code}\n")
        return 1
    except URLError as error:
        sys.stderr.write(f"Network error while fetching record {stable_id}: {error}\n")
        return 1

    write_json(output_dir / "brca1_record.json", full_record)

    try:
        genomic_fasta = client.fetch_sequence_fasta(stable_id, seq_type="genomic")
    except HTTPError as error:
        sys.stderr.write(f"Failed to fetch genomic sequence for {stable_id}: HTTP {error.code}\n")
        return 1
    except URLError as error:
        sys.stderr.write(f"Network error while fetching genomic sequence {stable_id}: {error}\n")
        return 1

    write_text(output_dir / "brca1_genomic.fasta", genomic_fasta)

    transcript_id = full_record.get("canonical_transcript")
    if transcript_id:
        try:
            cdna_fasta = client.fetch_sequence_fasta(transcript_id, seq_type="cdna")
            write_text(output_dir / "brca1_canonical_cdna.fasta", cdna_fasta)
        except HTTPError as error:
            sys.stderr.write(f"Failed to fetch cDNA sequence for {transcript_id}: HTTP {error.code}\n")
        except URLError as error:
            sys.stderr.write(f"Network error while fetching cDNA sequence {transcript_id}: {error}\n")

    protein_id = full_record.get("canonical_protein")
    if protein_id and not args.skip_protein:
        try:
            protein_fasta = client.fetch_sequence_fasta(protein_id, seq_type="protein")
            write_text(output_dir / "brca1_canonical_protein.fasta", protein_fasta)
        except HTTPError as error:
            sys.stderr.write(f"Failed to fetch protein sequence for {protein_id}: HTTP {error.code}\n")
        except URLError as error:
            sys.stderr.write(f"Network error while fetching protein sequence {protein_id}: {error}\n")

    print(f"Fetched BRCA1 data for {stable_id} into {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

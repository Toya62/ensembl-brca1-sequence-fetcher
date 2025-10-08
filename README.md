# Ensembl BRCA1 Sequence Fetcher

This repository contains a small Python utility that looks up the human BRCA1 gene in Ensembl, downloads the associated metadata, and saves canonical sequence data as FASTA files. It uses the Ensembl REST API (https://rest.ensembl.org/) and produces results in the `outputs/` directory by default.

## Usage

```bash
python3 fetch_brca1.py --output outputs
```

The script writes the following artifacts:
- `brca1_lookup.json`: raw response from the lookup-by-symbol endpoint.
- `brca1_record.json`: expanded record for the BRCA1 Ensembl stable ID.
- `brca1_genomic.fasta`: genomic DNA sequence for the gene.
- `brca1_canonical_cdna.fasta`: cDNA sequence for the canonical transcript (if available).
- `brca1_canonical_protein.fasta`: protein sequence for the canonical translation (can be skipped with `--skip-protein`).

Use `--symbol` and `--species` if you need to target a different gene or organism.

> **Note:** If you encounter TLS certificate issues in restricted environments, pass `--insecure` to disable verification temporarily.

#!/usr/bin/env python3
"""Parse strand-aware gene transcription start sites (TSS) from a GENCODE GTF file."""

from __future__ import annotations

import argparse
import gzip
import logging
import re
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse GTF key-value attribute string into a dictionary."""
    return dict(ATTR_RE.findall(attr_str))


def extract_tss_table(
    gtf_path: Path,
    chromosomes: set[str] | None = None,
    biotypes: set[str] | None = None,
) -> pd.DataFrame:
    """Extract gene-level strand-aware TSS coordinates from GENCODE GTF."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    records: list[dict[str, object]] = []

    with opener(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            chrom = parts[0]
            if chromosomes and chrom not in chromosomes:
                continue

            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parse_attributes(parts[8])

            biotype = attrs.get("gene_type") or attrs.get("gene_biotype", "")
            if biotypes and biotype not in biotypes:
                continue

            gene_id_raw = attrs.get("gene_id", "")
            gene_id = gene_id_raw.split(".")[0]
            gene_name = attrs.get("gene_name", gene_id)
            tss = start if strand == "+" else end

            records.append(
                {
                    "chrom": chrom,
                    "tss": tss,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_id_raw": gene_id_raw,
                    "gene_name": gene_name,
                    "gene_biotype": biotype,
                }
            )

    df = pd.DataFrame(records)
    if df.empty:
        raise ValueError(f"No gene TSS records extracted from {gtf_path}.")

    df = (
        df.sort_values(["chrom", "tss", "gene_id"])
        .drop_duplicates(subset=["gene_id"], keep="first")
        .reset_index(drop=True)
    )
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract strand-aware gene TSSs from GENCODE GTF.")
    parser.add_argument("--gtf", required=True, type=Path, help="Path to GENCODE .gtf or .gtf.gz.")
    parser.add_argument("--out", required=True, type=Path, help="Output TSV path.")
    parser.add_argument("--chromosomes", default="", help="Comma-separated chromosome whitelist.")
    parser.add_argument("--biotypes", default="protein_coding,lincRNA,antisense", help="Comma-separated gene biotypes.")
    args = parser.parse_args()

    chrom_set = {c.strip() for c in args.chromosomes.split(",") if c.strip()} or None
    bio_set = {b.strip() for b in args.biotypes.split(",") if b.strip()} or None

    df = extract_tss_table(args.gtf, chromosomes=chrom_set, biotypes=bio_set)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    logger.info("Extracted %d gene TSSs -> %s", len(df), args.out)


if __name__ == "__main__":
    main()

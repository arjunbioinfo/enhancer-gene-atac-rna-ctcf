#!/usr/bin/env python3
"""Download paired public ENCODE GRCh38 ATAC-seq, CTCF ChIP-seq, RNA-seq, and GENCODE GTF."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd
import requests
import yaml

ENCODE_BASE = "https://www.encodeproject.org"
HEADERS = {"Accept": "application/json", "User-Agent": "enhancer-gene-atac-rna-ctcf/1.0"}

PEAK_OUTPUT_PRIORITY = [
    "IDR thresholded peaks",
    "conservative IDR thresholded peaks",
    "optimal IDR thresholded peaks",
    "pseudoreplicated peaks",
    "replicated peaks",
    "peaks",
]

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def fetch_experiment_json(accession: str, timeout: int = 30) -> dict[str, Any]:
    """Query the ENCODE REST API for an experiment accession."""
    url = f"{ENCODE_BASE}/experiments/{accession}/?format=json"
    logger.info("Querying ENCODE API: %s", url)
    resp = requests.get(url, headers=HEADERS, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def select_peak_file(exp_json: dict[str, Any], assembly: str = "GRCh38") -> dict[str, Any]:
    """Select the highest-priority released narrowPeak file for the requested assembly."""
    candidates = [
        f
        for f in exp_json.get("files", [])
        if f.get("status") == "released"
        and f.get("assembly") == assembly
        and f.get("file_format") == "bed"
        and f.get("file_format_type") == "narrowPeak"
    ]
    if not candidates:
        raise ValueError(
            f"No released {assembly} bed narrowPeak files found for {exp_json.get('accession')}."
        )

    def sort_key(f: dict[str, Any]) -> tuple[int, int]:
        out_type = f.get("output_type", "")
        prio = (
            PEAK_OUTPUT_PRIORITY.index(out_type)
            if out_type in PEAK_OUTPUT_PRIORITY
            else len(PEAK_OUTPUT_PRIORITY)
        )
        # Prefer combined multi-replicate files (more biological replicates)
        n_reps = len(f.get("biological_replicates", []))
        return (prio, -n_reps)

    candidates.sort(key=sort_key)
    return candidates[0]


def select_rna_quant_files(
    exp_json: dict[str, Any], assembly: str = "GRCh38"
) -> list[dict[str, Any]]:
    """Select released gene quantification TSV files for the requested assembly."""
    candidates = [
        f
        for f in exp_json.get("files", [])
        if f.get("status") == "released"
        and f.get("assembly") == assembly
        and f.get("file_format") == "tsv"
        and f.get("output_type") == "gene quantifications"
    ]
    if not candidates:
        raise ValueError(
            f"No released {assembly} gene quantifications TSV files found for {exp_json.get('accession')}."
        )

    # Prefer V29 annotation if available so ENSG IDs align with GENCODE v29
    v29 = [f for f in candidates if f.get("genome_annotation") == "V29"]
    selected = v29 if v29 else candidates
    selected.sort(key=lambda f: f.get("accession", ""))
    return selected[:2]  # Up to 2 biological replicates for fast download


def download_file(url: str, dest: Path, timeout: int = 120) -> int:
    """Stream-download a remote URL to disk if not already present."""
    if dest.exists() and dest.stat().st_size > 0:
        logger.info("Using cached file: %s (%d bytes)", dest, dest.stat().st_size)
        return dest.stat().st_size

    dest.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %s -> %s", url, dest)
    tmp_path = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=timeout, headers={"User-Agent": HEADERS["User-Agent"]}) as r:
        r.raise_for_status()
        with open(tmp_path, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1024 * 256):
                if chunk:
                    fh.write(chunk)
    tmp_path.rename(dest)
    size = dest.stat().st_size
    logger.info("Saved %s (%d bytes)", dest, size)
    return size


def merge_rna_replicates(tsv_paths: list[Path], out_path: Path) -> None:
    """Average TPM and expected_count across ENCODE RSEM gene quantification replicates."""
    dfs = []
    for p in tsv_paths:
        df = pd.read_csv(p, sep="\t", usecols=lambda c: c in {"gene_id", "TPM", "expected_count", "FPKM"})
        if "gene_id" not in df.columns or "TPM" not in df.columns:
            raise ValueError(f"Unexpected RNA quantification columns in {p}: {list(df.columns)}")
        df["gene_id_clean"] = df["gene_id"].astype(str).str.split(".").str[0]
        dfs.append(df[["gene_id_clean", "TPM", "expected_count"]])

    combined = pd.concat(dfs, ignore_index=True)
    agg = (
        combined.groupby("gene_id_clean", as_index=False)[["TPM", "expected_count"]]
        .mean()
        .rename(columns={"gene_id_clean": "gene_id"})
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote replicate-averaged RNA table: %s (%d genes)", out_path, len(agg))


def run_download(config_path: Path, out_dir: Path) -> None:
    """Download all ENCODE assay files and GENCODE GTF specified in config."""
    with open(config_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    assembly = cfg["genome"]["assembly"]
    cond_a = cfg["comparison"]["condition_a"]
    cond_b = cfg["comparison"]["condition_b"]
    manifest: dict[str, Any] = {"assembly": assembly, "files": {}}

    for cond in (cond_a, cond_b):
        enc_cfg = cfg["encode"][cond]

        # 1. ATAC-seq narrowPeak
        atac_exp = fetch_experiment_json(enc_cfg["atac_experiment"])
        atac_file = select_peak_file(atac_exp, assembly=assembly)
        atac_url = ENCODE_BASE + atac_file["href"]
        atac_dest = out_dir / f"{cond}_atac.narrowPeak.gz"
        size = download_file(atac_url, atac_dest)
        manifest["files"][f"{cond}_atac"] = {
            "experiment": enc_cfg["atac_experiment"],
            "accession": atac_file["accession"],
            "output_type": atac_file.get("output_type"),
            "url": atac_url,
            "bytes": size,
        }

        # 2. CTCF ChIP-seq narrowPeak
        ctcf_exp = fetch_experiment_json(enc_cfg["ctcf_experiment"])
        ctcf_file = select_peak_file(ctcf_exp, assembly=assembly)
        ctcf_url = ENCODE_BASE + ctcf_file["href"]
        ctcf_dest = out_dir / f"{cond}_ctcf.narrowPeak.gz"
        size = download_file(ctcf_url, ctcf_dest)
        manifest["files"][f"{cond}_ctcf"] = {
            "experiment": enc_cfg["ctcf_experiment"],
            "accession": ctcf_file["accession"],
            "output_type": ctcf_file.get("output_type"),
            "url": ctcf_url,
            "bytes": size,
        }

        # 3. RNA-seq gene quantifications
        rna_exp = fetch_experiment_json(enc_cfg["rna_experiment"])
        rna_files = select_rna_quant_files(rna_exp, assembly=assembly)
        rep_paths: list[Path] = []
        rep_meta = []
        for idx, rf in enumerate(rna_files, start=1):
            r_url = ENCODE_BASE + rf["href"]
            r_dest = out_dir / f"{cond}_rna_rep{idx}_{rf['accession']}.tsv"
            r_size = download_file(r_url, r_dest)
            rep_paths.append(r_dest)
            rep_meta.append(
                {
                    "accession": rf["accession"],
                    "genome_annotation": rf.get("genome_annotation"),
                    "url": r_url,
                    "bytes": r_size,
                }
            )
        merge_rna_replicates(rep_paths, out_dir / f"{cond}_rna_quant.tsv")
        manifest["files"][f"{cond}_rna"] = {
            "experiment": enc_cfg["rna_experiment"],
            "replicates": rep_meta,
        }

    # 4. GENCODE basic GTF
    gtf_url = cfg["genome"]["gencode_gtf_url"]
    gtf_dest = out_dir / "gencode.basic.annotation.gtf.gz"
    gtf_size = download_file(gtf_url, gtf_dest)
    manifest["files"]["gencode_gtf"] = {"url": gtf_url, "bytes": gtf_size}

    manifest_path = out_dir / "download_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    logger.info("Wrote provenance manifest: %s", manifest_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Download public ENCODE and GENCODE data.")
    parser.add_argument("--config", required=True, type=Path, help="Path to YAML config file.")
    parser.add_argument("--out-dir", required=True, type=Path, help="Output directory for raw files.")
    args = parser.parse_args()
    run_download(args.config, args.out_dir)


if __name__ == "__main__":
    main()

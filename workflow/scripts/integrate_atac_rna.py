#!/usr/bin/env python3
"""Integrate condition-biased ATAC-seq accessibility with RNA-seq expression across gene TSSs.

Accessibility and expression changes are fold-change screens between one ENCODE experiment per
condition. They flag condition-biased loci; they are not replicate-based differential tests."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

NARROWPEAK_COLS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak",
]


def load_narrowpeak(path: Path, chromosomes: set[str] | None = None) -> pd.DataFrame:
    """Load an ENCODE narrowPeak (.bed or .bed.gz) file."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=range(10),
        names=NARROWPEAK_COLS,
        compression="infer",
    )
    if chromosomes:
        df = df[df["chrom"].isin(chromosomes)].copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["signalValue"] = pd.to_numeric(df["signalValue"], errors="coerce").fillna(0.0)
    df["peak"] = pd.to_numeric(df["peak"], errors="coerce").fillna(-1).astype(int)
    return df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)


def peak_summits(peaks: pd.DataFrame) -> np.ndarray:
    """Absolute summit position: start + narrowPeak offset (column 10), or the midpoint if absent (-1)."""
    starts = peaks["start"].to_numpy(dtype=int)
    ends = peaks["end"].to_numpy(dtype=int)
    if "peak" not in peaks.columns:
        return (starts + ends) // 2
    offsets = pd.to_numeric(peaks["peak"], errors="coerce").fillna(-1).to_numpy(dtype=int)
    return np.where(offsets >= 0, starts + offsets, (starts + ends) // 2)


def build_consensus_peaks(
    peaks_a: pd.DataFrame,
    peaks_b: pd.DataFrame,
    cond_a: str,
    cond_b: str,
    atac_lfc_thresh: float = 1.0,
) -> pd.DataFrame:
    """Merge overlapping peaks across two conditions and flag condition-biased accessibility.

    Each consensus region keeps the summit of its strongest contributing peak. The log2 ratio
    compares peak-caller signal values from one experiment per condition; it is a descriptive
    fold-change screen, not a replicate-based test of differential accessibility.
    """
    a_tagged = peaks_a[["chrom", "start", "end", "signalValue"]].copy()
    a_tagged["peak_summit"] = peak_summits(peaks_a)
    a_tagged["source"] = "A"
    b_tagged = peaks_b[["chrom", "start", "end", "signalValue"]].copy()
    b_tagged["peak_summit"] = peak_summits(peaks_b)
    b_tagged["source"] = "B"

    combined = pd.concat([a_tagged, b_tagged], ignore_index=True)
    combined = combined.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

    merged_rows: list[dict[str, object]] = []
    for chrom, grp in combined.groupby("chrom", sort=True):
        starts = grp["start"].to_numpy()
        ends = grp["end"].to_numpy()
        sigs = grp["signalValue"].to_numpy()
        srcs = grp["source"].to_numpy()
        pk_summits = grp["peak_summit"].to_numpy()

        def new_cluster(i: int) -> dict[str, object]:
            return {
                "start": int(starts[i]),
                "end": int(ends[i]),
                "max_a": float(sigs[i]) if srcs[i] == "A" else 0.0,
                "max_b": float(sigs[i]) if srcs[i] == "B" else 0.0,
                "in_a": srcs[i] == "A",
                "in_b": srcs[i] == "B",
                "best_sig": float(sigs[i]),
                "summit": int(pk_summits[i]),
            }

        def close_cluster(c: dict[str, object]) -> None:
            merged_rows.append(
                {
                    "chrom": chrom,
                    "start": c["start"],
                    "end": c["end"],
                    "summit": c["summit"],
                    f"atac_signal_{cond_a}": c["max_a"],
                    f"atac_signal_{cond_b}": c["max_b"],
                    f"present_{cond_a}": c["in_a"],
                    f"present_{cond_b}": c["in_b"],
                }
            )

        cur = new_cluster(0)
        for i in range(1, len(starts)):
            # BED intervals are half-open: a peak starting exactly at cur_end is adjacent, not overlapping
            if int(starts[i]) < cur["end"]:
                cur["end"] = max(cur["end"], int(ends[i]))
                sig = float(sigs[i])
                if srcs[i] == "A":
                    cur["in_a"] = True
                    cur["max_a"] = max(cur["max_a"], sig)
                else:
                    cur["in_b"] = True
                    cur["max_b"] = max(cur["max_b"], sig)
                if sig > cur["best_sig"]:
                    cur["best_sig"] = sig
                    cur["summit"] = int(pk_summits[i])
            else:
                close_cluster(cur)
                cur = new_cluster(i)
        close_cluster(cur)

    consensus = pd.DataFrame(merged_rows)
    consensus.insert(0, "peak_id", [f"PEAK_{i + 1:06d}" for i in range(len(consensus))])

    sig_a = consensus[f"atac_signal_{cond_a}"].to_numpy()
    sig_b = consensus[f"atac_signal_{cond_b}"].to_numpy()
    consensus["atac_log2fc"] = np.log2((sig_a + 1.0) / (sig_b + 1.0))

    conditions = [
        consensus["atac_log2fc"] >= atac_lfc_thresh,
        consensus["atac_log2fc"] <= -atac_lfc_thresh,
    ]
    choices = [f"{cond_a}_specific", f"{cond_b}_specific"]
    consensus["atac_class"] = np.select(conditions, choices, default="shared")
    return consensus


def link_peaks_to_genes(
    consensus: pd.DataFrame,
    tss_df: pd.DataFrame,
    promoter_bp: int = 2000,
    enhancer_max_bp: int = 250000,
    decay_scale_bp: float = 50000.0,
) -> pd.DataFrame:
    """Link consensus ATAC peaks to gene TSSs within enhancer_max_bp using binary search."""
    link_rows: list[dict[str, object]] = []

    tss_by_chrom = {
        chrom: grp.sort_values("tss").reset_index(drop=True)
        for chrom, grp in tss_df.groupby("chrom")
    }

    for chrom, p_grp in consensus.groupby("chrom", sort=True):
        if chrom not in tss_by_chrom:
            continue
        g_grp = tss_by_chrom[chrom]
        tss_coords = g_grp["tss"].to_numpy()
        gene_ids = g_grp["gene_id"].to_numpy()
        gene_names = g_grp["gene_name"].to_numpy()
        strands = g_grp["strand"].to_numpy()

        for row in p_grp.itertuples(index=False):
            summit = int(row.summit)
            left = int(np.searchsorted(tss_coords, summit - enhancer_max_bp, side="left"))
            right = int(np.searchsorted(tss_coords, summit + enhancer_max_bp, side="right"))
            if left >= right:
                continue

            cand_tss = tss_coords[left:right]
            signed_dists = summit - cand_tss
            abs_dists = np.abs(signed_dists)
            nearest_local_idx = int(np.argmin(abs_dists))

            for offset in range(right - left):
                g_idx = left + offset
                abs_d = int(abs_dists[offset])
                signed_d = int(signed_dists[offset])
                if strands[g_idx] == "-":
                    signed_d = -signed_d

                elem_type = "promoter" if abs_d <= promoter_bp else "distal_enhancer"
                weight = float(np.exp(-abs_d / decay_scale_bp))

                link_rows.append(
                    {
                        "peak_id": row.peak_id,
                        "chrom": chrom,
                        "start": int(row.start),
                        "end": int(row.end),
                        "summit": summit,
                        "atac_log2fc": float(row.atac_log2fc),
                        "atac_class": row.atac_class,
                        "gene_id": gene_ids[g_idx],
                        "gene_name": gene_names[g_idx],
                        "tss": int(tss_coords[g_idx]),
                        "signed_dist_bp": signed_d,
                        "abs_dist_bp": abs_d,
                        "element_type": elem_type,
                        "is_nearest_gene": offset == nearest_local_idx,
                        "distance_weight": weight,
                        "weighted_atac_log2fc": float(row.atac_log2fc) * weight,
                    }
                )

    return pd.DataFrame(link_rows)


def integrate_expression(
    links: pd.DataFrame,
    tss_df: pd.DataFrame,
    rna_a: pd.DataFrame,
    rna_b: pd.DataFrame,
    cond_a: str,
    cond_b: str,
    min_tpm: float = 1.0,
    deg_lfc_thresh: float = 1.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join RNA-seq expression onto enhancer-gene links and build per-gene summary."""
    expr = tss_df[["gene_id", "gene_name", "chrom", "tss", "gene_biotype"]].merge(
        rna_a[["gene_id", "TPM"]].rename(columns={"TPM": f"tpm_{cond_a}"}),
        on="gene_id",
        how="left",
    ).merge(
        rna_b[["gene_id", "TPM"]].rename(columns={"TPM": f"tpm_{cond_b}"}),
        on="gene_id",
        how="left",
    )
    expr[f"tpm_{cond_a}"] = expr[f"tpm_{cond_a}"].fillna(0.0)
    expr[f"tpm_{cond_b}"] = expr[f"tpm_{cond_b}"].fillna(0.0)
    expr["expressed"] = (expr[f"tpm_{cond_a}"] >= min_tpm) | (expr[f"tpm_{cond_b}"] >= min_tpm)
    expr["rna_log2fc"] = np.log2((expr[f"tpm_{cond_a}"] + 1.0) / (expr[f"tpm_{cond_b}"] + 1.0))

    deg_conds = [
        expr["expressed"] & (expr["rna_log2fc"] >= deg_lfc_thresh),
        expr["expressed"] & (expr["rna_log2fc"] <= -deg_lfc_thresh),
    ]
    expr["deg_status"] = np.select(deg_conds, [f"{cond_a}_up", f"{cond_b}_up"], default="unchanged")

    links_merged = links.merge(
        expr[["gene_id", f"tpm_{cond_a}", f"tpm_{cond_b}", "expressed", "rna_log2fc", "deg_status"]],
        on="gene_id",
        how="left",
    )

    # Aggregate enhancer and promoter signals per gene
    distal = links_merged[links_merged["element_type"] == "distal_enhancer"]
    prom = links_merged[links_merged["element_type"] == "promoter"]

    distal_agg = (
        distal.groupby("gene_id", as_index=False)
        .agg(
            n_distal_enhancers=("peak_id", "nunique"),
            n_enhancers_a_spec=("atac_class", lambda s: int((s == f"{cond_a}_specific").sum())),
            n_enhancers_b_spec=("atac_class", lambda s: int((s == f"{cond_b}_specific").sum())),
            mean_enhancer_atac_lfc=("atac_log2fc", "mean"),
            weighted_enhancer_score=("weighted_atac_log2fc", "sum"),
        )
    )
    prom_agg = (
        prom.groupby("gene_id", as_index=False)
        .agg(
            n_promoter_peaks=("peak_id", "nunique"),
            mean_promoter_atac_lfc=("atac_log2fc", "mean"),
        )
    )

    gene_summary = (
        expr.merge(distal_agg, on="gene_id", how="left")
        .merge(prom_agg, on="gene_id", how="left")
    )
    for col in ("n_distal_enhancers", "n_enhancers_a_spec", "n_enhancers_b_spec", "n_promoter_peaks"):
        gene_summary[col] = gene_summary[col].fillna(0).astype(int)
    for col in ("mean_enhancer_atac_lfc", "weighted_enhancer_score", "mean_promoter_atac_lfc"):
        gene_summary[col] = gene_summary[col].fillna(0.0)

    # Classify multi-omic regulatory concordance
    conc_conds = [
        (gene_summary["deg_status"] == f"{cond_a}_up") & (gene_summary["weighted_enhancer_score"] > 0.5),
        (gene_summary["deg_status"] == f"{cond_b}_up") & (gene_summary["weighted_enhancer_score"] < -0.5),
        (gene_summary["deg_status"] == "unchanged") & (gene_summary["weighted_enhancer_score"].abs() >= 1.5),
        (
            (gene_summary["deg_status"] == f"{cond_a}_up") & (gene_summary["weighted_enhancer_score"] < -0.5)
        )
        | (
            (gene_summary["deg_status"] == f"{cond_b}_up") & (gene_summary["weighted_enhancer_score"] > 0.5)
        ),
    ]
    gene_summary["regulatory_class"] = np.select(
        conc_conds,
        [f"concordant_{cond_a}", f"concordant_{cond_b}", "chromatin_primed", "discordant"],
        default="baseline",
    )

    return links_merged, gene_summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate ATAC-seq peaks with RNA-seq expression.")
    parser.add_argument("--atac-a", required=True, type=Path)
    parser.add_argument("--atac-b", required=True, type=Path)
    parser.add_argument("--rna-a", required=True, type=Path)
    parser.add_argument("--rna-b", required=True, type=Path)
    parser.add_argument("--tss", required=True, type=Path)
    parser.add_argument("--cond-a", required=True)
    parser.add_argument("--cond-b", required=True)
    parser.add_argument("--chromosomes", default="")
    parser.add_argument("--promoter-bp", type=int, default=2000)
    parser.add_argument("--enhancer-max-bp", type=int, default=250000)
    parser.add_argument("--min-tpm", type=float, default=1.0)
    parser.add_argument("--deg-lfc", type=float, default=1.0)
    parser.add_argument("--atac-lfc", type=float, default=1.0)
    parser.add_argument("--out-links", required=True, type=Path)
    parser.add_argument("--out-genes", required=True, type=Path)
    args = parser.parse_args()

    chrom_set = {c.strip() for c in args.chromosomes.split(",") if c.strip()} or None

    peaks_a = load_narrowpeak(args.atac_a, chromosomes=chrom_set)
    peaks_b = load_narrowpeak(args.atac_b, chromosomes=chrom_set)
    tss_df = pd.read_csv(args.tss, sep="\t")
    rna_a = pd.read_csv(args.rna_a, sep="\t")
    rna_b = pd.read_csv(args.rna_b, sep="\t")

    consensus = build_consensus_peaks(peaks_a, peaks_b, args.cond_a, args.cond_b, args.atac_lfc)
    logger.info("Built %d consensus ATAC peaks.", len(consensus))

    links = link_peaks_to_genes(consensus, tss_df, args.promoter_bp, args.enhancer_max_bp)
    links_merged, gene_summary = integrate_expression(
        links, tss_df, rna_a, rna_b, args.cond_a, args.cond_b, args.min_tpm, args.deg_lfc
    )

    args.out_links.parent.mkdir(parents=True, exist_ok=True)
    links_merged.to_csv(args.out_links, sep="\t", index=False)
    gene_summary.to_csv(args.out_genes, sep="\t", index=False)
    logger.info("Saved %d peak-gene links and %d gene summaries.", len(links_merged), len(gene_summary))


if __name__ == "__main__":
    main()

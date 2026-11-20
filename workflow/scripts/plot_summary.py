#!/usr/bin/env python3
"""Generate publication-quality summary figures for ATAC-seq, RNA-seq, and CTCF integration."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)


def plot_atac_rna_correlation(
    genes: pd.DataFrame, cond_a: str, cond_b: str, highlight_genes: list[str], out_path: Path
) -> None:
    """Scatter plot of distance-weighted distal enhancer Δaccessibility vs RNA-seq log2FC."""
    df = genes[genes["expressed"] & (genes["n_distal_enhancers"] > 0)].copy()
    fig, ax = plt.subplots(figsize=(7.2, 5.6))

    palette = {
        f"concordant_{cond_a}": "#D95F02",
        f"concordant_{cond_b}": "#1B9E77",
        "chromatin_primed": "#7570B3",
        "discordant": "#E7298A",
        "baseline": "#BDBDBD",
    }

    for cls, grp in df.groupby("regulatory_class"):
        ax.scatter(
            grp["weighted_enhancer_score"],
            grp["rna_log2fc"],
            s=18 if cls == "baseline" else 30,
            alpha=0.45 if cls == "baseline" else 0.8,
            color=palette.get(str(cls), "#999999"),
            label=f"{cls} (n={len(grp)})",
            edgecolors="none",
        )

    ax.axhline(0, color="#444444", lw=0.8, ls="--")
    ax.axvline(0, color="#444444", lw=0.8, ls="--")

    if len(df) >= 3:
        rho, pval = spearmanr(df["weighted_enhancer_score"], df["rna_log2fc"])
        ax.text(
            0.03,
            0.96,
            f"Spearman ρ = {rho:.3f} (p = {pval:.2e})\nn = {len(df)} expressed genes",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"},
        )

    # Annotate requested highlight genes if present
    hl = df[df["gene_name"].isin(highlight_genes)]
    for row in hl.itertuples(index=False):
        ax.scatter(
            row.weighted_enhancer_score,
            row.rna_log2fc,
            s=65,
            color="black",
            zorder=5,
            facecolors="none",
            linewidths=1.4,
        )
        ax.annotate(
            row.gene_name,
            xy=(row.weighted_enhancer_score, row.rna_log2fc),
            xytext=(6, 5),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
            color="#111111",
        )

    ax.set_xlabel(f"Distance-weighted distal enhancer ΔATAC ({cond_a} / {cond_b})")
    ax.set_ylabel(f"RNA-seq log2FC ({cond_a} / {cond_b})")
    ax.set_title(f"Distal Enhancer Accessibility vs. Gene Expression ({cond_a} vs. {cond_b})")
    ax.legend(frameon=True, fontsize=8.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_ctcf_enrichment(stats_df: pd.DataFrame, out_path: Path) -> None:
    """Forest plot of CTCF occupancy odds ratios with 95% confidence intervals."""
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    df = stats_df.copy()
    if "test" in df.columns:
        df = df[df["test"] == "fisher"]
    df = df.reset_index(drop=True)
    y_pos = np.arange(len(df))

    ors = np.clip(df["odds_ratio"].to_numpy(dtype=float), 0.05, 50.0)
    ci_l = np.clip(df["ci_low"].to_numpy(dtype=float), 0.05, 50.0)
    ci_h = np.clip(df["ci_high"].to_numpy(dtype=float), 0.05, 50.0)
    xerr = np.vstack([np.maximum(ors - ci_l, 0.0), np.maximum(ci_h - ors, 0.0)])

    colors = ["#D95F02" if fdr < 0.05 else "#636363" for fdr in df["fdr_bh"]]
    ax.errorbar(ors, y_pos, xerr=xerr, fmt="o", color="#222222", ecolor="#555555", capsize=4, lw=1.4, zorder=3)
    ax.scatter(ors, y_pos, c=colors, s=70, zorder=4)
    ax.axvline(1.0, color="#999999", ls="--", lw=1.0)
    ax.set_xscale("log")

    labels = [c.replace("_", " ") for c in df["comparison"]]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    for idx, row in df.iterrows():
        ax.text(
            ci_h[idx] * 1.08,
            idx,
            f"OR={row['odds_ratio']:.2f}, FDR={row['fdr_bh']:.2e}",
            va="center",
            fontsize=8.5,
        )

    if len(df):
        lo, hi = float(np.min(ci_l)), float(np.max(ci_h))
        ax.set_xlim(min(lo, 1.0) / 1.4, max(hi, 1.0) * 4.0)
        ticks = [t for t in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50) if ax.get_xlim()[0] <= t <= ax.get_xlim()[1]]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.minorticks_off()
    ax.set_xlabel("CTCF Occupancy Odds Ratio (log scale, 95% CI)")
    ax.set_title("CTCF Occupancy at Condition-Biased Enhancers (peak-level)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_distance_to_ctcf(links: pd.DataFrame, cond_a: str, cond_b: str, out_path: Path) -> None:
    """Empirical CDF of distance from distal enhancers to nearest CTCF peak."""
    distal = links[links["element_type"] == "distal_enhancer"].drop_duplicates(subset=["peak_id"]).copy()
    fig, ax = plt.subplots(figsize=(6.8, 4.8))

    palette = {
        f"{cond_a}_specific": "#D95F02",
        f"{cond_b}_specific": "#1B9E77",
        "shared": "#386CB0",
    }

    for cls in (f"{cond_a}_specific", f"{cond_b}_specific", "shared"):
        vals = distal[distal["atac_class"] == cls]["dist_ctcf_min_bp"].dropna().to_numpy()
        if len(vals) == 0:
            continue
        sorted_kb = np.sort(vals) / 1000.0
        ecdf = np.arange(1, len(sorted_kb) + 1) / len(sorted_kb)
        ax.plot(sorted_kb, ecdf, lw=2.0, color=palette.get(cls, "#555555"), label=f"{cls} (n={len(vals)})")

    ax.set_xscale("log")
    ax.set_xlim(0.05, 500)
    ax.set_xlabel("Distance to nearest CTCF ChIP-seq summit (kb, log scale)")
    ax.set_ylabel("Cumulative fraction of distal enhancers")
    ax.set_title("Proximity of Distal Enhancers to CTCF Binding Sites")
    ax.legend(frameon=True, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_top_candidates(candidates: pd.DataFrame, cond_a: str, cond_b: str, out_path: Path) -> None:
    """Horizontal bar/dot summary of top 15 prioritized enhancer-gene-CTCF candidates."""
    top = candidates.head(15).copy().iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 5.4))

    if top.empty:
        ax.text(0.5, 0.5, "No candidate genes passed filters.", ha="center", va="center")
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        return

    colors = ["#D95F02" if lfc > 0 else "#1B9E77" for lfc in top["rna_log2fc"]]
    ax.barh(top["gene_name"], top["priority_score"], color=colors, alpha=0.85, height=0.65)

    for idx, row in enumerate(top.itertuples(index=False)):
        ax.text(
            row.priority_score + max(top["priority_score"]) * 0.015,
            idx,
            f"RNA LFC={row.rna_log2fc:+.1f} | {row.n_ctcf_bound_enhancers}/{row.n_distal_enhancers} CTCF+ enh",
            va="center",
            fontsize=8.2,
        )

    ax.set_xlabel("Multi-Omic Regulatory Priority Score")
    ax.set_title(f"Top Prioritized Enhancer–Gene–CTCF Candidates ({cond_a} orange / {cond_b} teal)")
    ax.set_xlim(0, max(top["priority_score"]) * 1.48)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot multi-omic ATAC/RNA/CTCF integration results.")
    parser.add_argument("--links", required=True, type=Path)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--cond-a", required=True)
    parser.add_argument("--cond-b", required=True)
    parser.add_argument("--highlight-genes", default="MYC,PVT1,GATA1,CD79A,PAX5")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    links = pd.read_csv(args.links, sep="\t")
    genes = pd.read_csv(args.genes, sep="\t")
    stats_df = pd.read_csv(args.stats, sep="\t")
    candidates = pd.read_csv(args.candidates, sep="\t")
    hl = [g.strip() for g in args.highlight_genes.split(",") if g.strip()]

    plot_atac_rna_correlation(genes, args.cond_a, args.cond_b, hl, args.out_dir / "atac_rna_correlation.png")
    plot_ctcf_enrichment(stats_df, args.out_dir / "ctcf_enrichment.png")
    plot_distance_to_ctcf(links, args.cond_a, args.cond_b, args.out_dir / "distance_to_ctcf.png")
    plot_top_candidates(candidates, args.cond_a, args.cond_b, args.out_dir / "top_candidates_summary.png")
    logger.info("Saved 4 summary figures to %s", args.out_dir)


if __name__ == "__main__":
    main()

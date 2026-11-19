#!/usr/bin/env python3
"""Annotate enhancer-gene links with CTCF ChIP-seq occupancy, test enrichment, and prioritize candidates."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# JASPAR MA0139.1 canonical 19-bp CTCF position frequency matrix (A, C, G, T)
# Used for optional sequence-level PWM scanning when a reference FASTA is provided.
CTCF_PFM_MA0139 = {
    "A": [87, 167, 281, 56, 8, 744, 40, 107, 851, 5, 333, 54, 12, 56, 104, 372, 82, 117, 402],
    "C": [291, 145, 49, 800, 903, 13, 528, 433, 11, 0, 3, 12, 528, 733, 13, 482, 322, 181, 162],
    "G": [76, 414, 449, 21, 0, 65, 334, 48, 32, 903, 566, 504, 290, 92, 750, 29, 437, 466, 200],
    "T": [459, 187, 134, 36, 2, 91, 11, 325, 19, 5, 11, 343, 83, 32, 46, 30, 72, 149, 149],
}


def score_ctcf_pwm(sequence: str) -> float:
    """Return the maximum log-odds score for the JASPAR MA0139.1 CTCF motif on both strands."""
    seq = sequence.upper()
    motif_len = len(CTCF_PFM_MA0139["A"])
    if len(seq) < motif_len:
        return float("-inf")

    bases = ("A", "C", "G", "T")
    comp = str.maketrans("ACGTN", "TGCAN")
    pfm = np.array([CTCF_PFM_MA0139[b] for b in bases], dtype=float)
    pwm = np.log2(((pfm + 0.25) / (pfm.sum(axis=0, keepdims=True) + 1.0)) / 0.25)
    base_to_idx = {b: i for i, b in enumerate(bases)}

    best_score = float("-inf")
    for strand_seq in (seq, seq.translate(comp)[::-1]):
        for i in range(len(strand_seq) - motif_len + 1):
            window = strand_seq[i : i + motif_len]
            if any(ch not in base_to_idx for ch in window):
                continue
            score = sum(pwm[base_to_idx[ch], pos] for pos, ch in enumerate(window))
            if score > best_score:
                best_score = float(score)
    return best_score


def load_ctcf_summits(path: Path, chromosomes: set[str] | None = None) -> dict[str, np.ndarray]:
    """Load CTCF narrowPeak file and return sorted summit coordinates per chromosome."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        compression="infer",
    )
    if chromosomes:
        df = df[df["chrom"].isin(chromosomes)].copy()
    df["summit"] = ((df["start"].astype(int) + df["end"].astype(int)) // 2).astype(int)
    return {
        chrom: np.sort(grp["summit"].to_numpy())
        for chrom, grp in df.groupby("chrom")
    }


def nearest_peak_distance(
    chroms: pd.Series, summits: pd.Series, target_by_chrom: dict[str, np.ndarray]
) -> np.ndarray:
    """Compute exact bp distance from each query summit to the nearest target summit."""
    dists = np.full(len(summits), np.nan, dtype=float)
    df_tmp = pd.DataFrame({"chrom": chroms.values, "summit": summits.values})

    for chrom, grp in df_tmp.groupby("chrom"):
        targets = target_by_chrom.get(str(chrom))
        if targets is None or len(targets) == 0:
            continue
        q = grp["summit"].to_numpy(dtype=int)
        idx = np.searchsorted(targets, q)
        d_right = np.where(idx < len(targets), np.abs(targets[np.minimum(idx, len(targets) - 1)] - q), np.inf)
        d_left = np.where(idx > 0, np.abs(q - targets[np.maximum(idx - 1, 0)]), np.inf)
        dists[grp.index.to_numpy()] = np.minimum(d_left, d_right)

    return dists


def run_fisher_test(fg_bound: int, fg_total: int, bg_bound: int, bg_total: int) -> dict[str, float]:
    """Two-sided Fisher's exact test with Haldane-Anscombe 95% confidence interval."""
    a = fg_bound
    b = max(fg_total - fg_bound, 0)
    c = bg_bound
    d = max(bg_total - bg_bound, 0)
    odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="two-sided")

    # Haldane-Anscombe continuity correction for log-OR SE
    a_c, b_c, c_c, d_c = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_or = np.log((a_c * d_c) / (b_c * c_c))
    se = np.sqrt(1.0 / a_c + 1.0 / b_c + 1.0 / c_c + 1.0 / d_c)
    ci_low = float(np.exp(log_or - 1.96 * se))
    ci_high = float(np.exp(log_or + 1.96 * se))

    return {
        "fg_bound": a,
        "fg_total": fg_total,
        "fg_frac": a / fg_total if fg_total else 0.0,
        "bg_bound": c,
        "bg_total": bg_total,
        "bg_frac": c / bg_total if bg_total else 0.0,
        "odds_ratio": float(odds_ratio),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "pvalue": float(pvalue),
    }


def annotate_and_test_ctcf(
    links: pd.DataFrame,
    genes: pd.DataFrame,
    ctcf_a: dict[str, np.ndarray],
    ctcf_b: dict[str, np.ndarray],
    cond_a: str,
    cond_b: str,
    ctcf_window_bp: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Annotate links with CTCF proximity, compute enrichment stats, and rank candidate genes."""
    out_links = links.copy()
    out_links[f"dist_ctcf_{cond_a}_bp"] = nearest_peak_distance(out_links["chrom"], out_links["summit"], ctcf_a)
    out_links[f"dist_ctcf_{cond_b}_bp"] = nearest_peak_distance(out_links["chrom"], out_links["summit"], ctcf_b)
    out_links["dist_ctcf_min_bp"] = np.nanmin(
        np.vstack([out_links[f"dist_ctcf_{cond_a}_bp"], out_links[f"dist_ctcf_{cond_b}_bp"]]),
        axis=0,
    )

    out_links[f"ctcf_bound_{cond_a}"] = out_links[f"dist_ctcf_{cond_a}_bp"] <= ctcf_window_bp
    out_links[f"ctcf_bound_{cond_b}"] = out_links[f"dist_ctcf_{cond_b}_bp"] <= ctcf_window_bp
    out_links["ctcf_bound_any"] = out_links[f"ctcf_bound_{cond_a}"] | out_links[f"ctcf_bound_{cond_b}"]

    # Deduplicate to unique distal enhancer peaks for peak-level enrichment tests
    distal_peaks = (
        out_links[out_links["element_type"] == "distal_enhancer"]
        .drop_duplicates(subset=["peak_id"])
        .copy()
    )

    tests: list[dict[str, object]] = []
    shared_bg = distal_peaks[distal_peaks["atac_class"] == "shared"]

    for label, subset, ctcf_col in [
        (f"{cond_a}_specific_enhancers_vs_shared", distal_peaks[distal_peaks["atac_class"] == f"{cond_a}_specific"], f"ctcf_bound_{cond_a}"),
        (f"{cond_b}_specific_enhancers_vs_shared", distal_peaks[distal_peaks["atac_class"] == f"{cond_b}_specific"], f"ctcf_bound_{cond_b}"),
    ]:
        res = run_fisher_test(
            int(subset[ctcf_col].sum()),
            len(subset),
            int(shared_bg[ctcf_col].sum()),
            len(shared_bg),
        )
        res["comparison"] = label
        tests.append(res)

    # Link-level enrichment: distal enhancers linked to DEGs vs unchanged genes
    distal_links = out_links[out_links["element_type"] == "distal_enhancer"]
    bg_deg = distal_links[distal_links["deg_status"] == "unchanged"]
    for label, deg_val, ctcf_col in [
        (f"enhancers_of_{cond_a}_up_DEGs_vs_unchanged", f"{cond_a}_up", f"ctcf_bound_{cond_a}"),
        (f"enhancers_of_{cond_b}_up_DEGs_vs_unchanged", f"{cond_b}_up", f"ctcf_bound_{cond_b}"),
    ]:
        fg = distal_links[distal_links["deg_status"] == deg_val]
        res = run_fisher_test(
            int(fg[ctcf_col].sum()),
            len(fg),
            int(bg_deg[ctcf_col].sum()),
            len(bg_deg),
        )
        res["comparison"] = label
        tests.append(res)

    stats_df = pd.DataFrame(tests)
    if not stats_df.empty:
        _, fdr, _, _ = multipletests(stats_df["pvalue"].to_numpy(), method="fdr_bh")
        stats_df["fdr_bh"] = fdr

        # Add Mann-Whitney U test for CTCF distance shift on condition-specific vs shared enhancers
        spec_d = distal_peaks[distal_peaks["atac_class"] != "shared"]["dist_ctcf_min_bp"].dropna()
        shar_d = shared_bg["dist_ctcf_min_bp"].dropna()
        if len(spec_d) > 0 and len(shar_d) > 0:
            mw_stat, mw_p = mannwhitneyu(spec_d, shar_d, alternative="two-sided")
            stats_df["mw_distance_pvalue"] = float(mw_p)
            stats_df["median_dist_specific_bp"] = float(np.median(spec_d))
            stats_df["median_dist_shared_bp"] = float(np.median(shar_d))

    # Build prioritized candidate table: genes with concordant enhancer + RNA shift + CTCF support
    ctcf_gene_agg = (
        distal_links.groupby("gene_id", as_index=False)
        .agg(
            n_ctcf_bound_enhancers=("ctcf_bound_any", "sum"),
            min_enhancer_ctcf_dist_bp=("dist_ctcf_min_bp", "min"),
        )
    )
    candidates = genes.merge(ctcf_gene_agg, on="gene_id", how="left")
    candidates["n_ctcf_bound_enhancers"] = candidates["n_ctcf_bound_enhancers"].fillna(0).astype(int)
    candidates["ctcf_enhancer_frac"] = np.where(
        candidates["n_distal_enhancers"] > 0,
        candidates["n_ctcf_bound_enhancers"] / candidates["n_distal_enhancers"],
        0.0,
    )
    # Composite priority score combining concordant RNA fold-change, weighted enhancer score, and CTCF support
    candidates["priority_score"] = (
        candidates["rna_log2fc"].abs()
        * np.log1p(candidates["weighted_enhancer_score"].abs())
        * (1.0 + 0.5 * candidates["ctcf_enhancer_frac"])
    )
    candidates = (
        candidates[candidates["expressed"] & (candidates["n_distal_enhancers"] > 0)]
        .sort_values("priority_score", ascending=False)
        .reset_index(drop=True)
    )

    return out_links, stats_df, candidates


def main() -> None:
    parser = argparse.ArgumentParser(description="CTCF ChIP-seq proximity and enrichment analysis.")
    parser.add_argument("--links", required=True, type=Path)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--ctcf-a", required=True, type=Path)
    parser.add_argument("--ctcf-b", required=True, type=Path)
    parser.add_argument("--cond-a", required=True)
    parser.add_argument("--cond-b", required=True)
    parser.add_argument("--chromosomes", default="")
    parser.add_argument("--ctcf-window-bp", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-stats", required=True, type=Path)
    parser.add_argument("--out-links", required=True, type=Path)
    parser.add_argument("--out-candidates", required=True, type=Path)
    args = parser.parse_args()

    np.random.seed(args.seed)
    chrom_set = {c.strip() for c in args.chromosomes.split(",") if c.strip()} or None

    links = pd.read_csv(args.links, sep="\t")
    genes = pd.read_csv(args.genes, sep="\t")
    ctcf_a = load_ctcf_summits(args.ctcf_a, chromosomes=chrom_set)
    ctcf_b = load_ctcf_summits(args.ctcf_b, chromosomes=chrom_set)

    out_links, stats_df, candidates = annotate_and_test_ctcf(
        links, genes, ctcf_a, ctcf_b, args.cond_a, args.cond_b, args.ctcf_window_bp
    )

    for p in (args.out_stats, args.out_links, args.out_candidates):
        p.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(args.out_stats, sep="\t", index=False)
    out_links.to_csv(args.out_links, sep="\t", index=False)
    candidates.to_csv(args.out_candidates, sep="\t", index=False)
    logger.info("Wrote CTCF enrichment stats (%d tests) and %d prioritized candidates.", len(stats_df), len(candidates))


if __name__ == "__main__":
    main()

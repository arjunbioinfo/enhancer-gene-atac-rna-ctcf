"""Unit and end-to-end smoke tests for the enhancer-gene ATAC/RNA/CTCF pipeline."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from ctcf_enrichment import annotate_and_test_ctcf, nearest_peak_distance, score_ctcf_pwm  # noqa: E402
from integrate_atac_rna import (  # noqa: E402
    build_consensus_peaks,
    integrate_expression,
    link_peaks_to_genes,
)
from parse_gencode_tss import extract_tss_table  # noqa: E402
from plot_summary import (  # noqa: E402
    plot_atac_rna_correlation,
    plot_ctcf_enrichment,
    plot_distance_to_ctcf,
    plot_top_candidates,
)


def test_parse_gencode_tss(tmp_path: Path) -> None:
    gtf_path = tmp_path / "mini.gtf.gz"
    lines = [
        "# comment line\n",
        'chr8\tHAVANA\tgene\t127735434\t127742951\t.\t+\t.\tgene_id "ENSG00000136997.16"; gene_type "protein_coding"; gene_name "MYC";\n',
        'chr8\tHAVANA\tgene\t127794533\t128101253\t.\t-\t.\tgene_id "ENSG00000249859.10"; gene_type "lincRNA"; gene_name "PVT1";\n',
        'chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000000001.1"; gene_type "pseudogene"; gene_name "SKIPME";\n',
    ]
    with gzip.open(gtf_path, "wt", encoding="utf-8") as fh:
        fh.writelines(lines)

    df = extract_tss_table(gtf_path, chromosomes={"chr8"}, biotypes={"protein_coding", "lincRNA"})
    assert list(df["gene_name"]) == ["MYC", "PVT1"]
    assert df.loc[df["gene_name"] == "MYC", "tss"].iloc[0] == 127735434
    # Minus-strand gene uses end coordinate as TSS
    assert df.loc[df["gene_name"] == "PVT1", "tss"].iloc[0] == 128101253


def test_score_ctcf_pwm_distinguishes_consensus() -> None:
    # Exact 19-bp JASPAR MA0139.1 CTCF consensus sequence
    canonical = "TGGCCACCAGGGCCGCGGA"
    poly_a = "AAAAAAAAAAAAAAAAAAA"
    assert score_ctcf_pwm(canonical) > 15.0
    assert score_ctcf_pwm(canonical) > score_ctcf_pwm(poly_a) + 15.0


def test_end_to_end_synthetic_integration_and_plots(tmp_path: Path) -> None:
    tss_df = pd.DataFrame(
        [
            {"chrom": "chr8", "tss": 100_000, "strand": "+", "gene_id": "ENSG_MYC", "gene_name": "MYC", "gene_biotype": "protein_coding"},
            {"chrom": "chr8", "tss": 400_000, "strand": "+", "gene_id": "ENSG_BCELL", "gene_name": "BLK", "gene_biotype": "protein_coding"},
            {"chrom": "chr8", "tss": 800_000, "strand": "-", "gene_id": "ENSG_CONST", "gene_name": "HOUSEKEEP", "gene_biotype": "protein_coding"},
        ]
    )

    peaks_a = pd.DataFrame(
        [
            {"chrom": "chr8", "start": 99_800, "end": 100_300, "signalValue": 40.0},   # MYC promoter
            {"chrom": "chr8", "start": 130_000, "end": 130_600, "signalValue": 55.0},  # MYC distal enhancer (K562-spec)
            {"chrom": "chr8", "start": 150_000, "end": 150_500, "signalValue": 35.0},  # MYC distal enhancer 2 (K562-spec)
            {"chrom": "chr8", "start": 790_000, "end": 790_500, "signalValue": 20.0},  # Shared distal enhancer
        ]
    )
    peaks_b = pd.DataFrame(
        [
            {"chrom": "chr8", "start": 99_900, "end": 100_250, "signalValue": 10.0},   # MYC promoter weaker
            {"chrom": "chr8", "start": 425_000, "end": 425_600, "signalValue": 60.0},  # BLK distal enhancer (GM-spec)
            {"chrom": "chr8", "start": 790_050, "end": 790_450, "signalValue": 21.0},  # Shared distal enhancer
        ]
    )

    consensus = build_consensus_peaks(peaks_a, peaks_b, "K562", "GM12878", atac_lfc_thresh=1.0)
    assert len(consensus) == 5
    assert (consensus["atac_class"] == "K562_specific").sum() >= 2
    assert (consensus["atac_class"] == "GM12878_specific").sum() == 1
    assert (consensus["atac_class"] == "shared").sum() == 1

    links = link_peaks_to_genes(consensus, tss_df, promoter_bp=2000, enhancer_max_bp=250_000)
    rna_a = pd.DataFrame(
        [
            {"gene_id": "ENSG_MYC", "TPM": 120.0},
            {"gene_id": "ENSG_BCELL", "TPM": 2.0},
            {"gene_id": "ENSG_CONST", "TPM": 30.0},
        ]
    )
    rna_b = pd.DataFrame(
        [
            {"gene_id": "ENSG_MYC", "TPM": 10.0},
            {"gene_id": "ENSG_BCELL", "TPM": 80.0},
            {"gene_id": "ENSG_CONST", "TPM": 28.0},
        ]
    )

    links_merged, gene_summary = integrate_expression(
        links, tss_df, rna_a, rna_b, "K562", "GM12878", min_tpm=1.0, deg_lfc_thresh=1.0
    )
    myc_row = gene_summary[gene_summary["gene_name"] == "MYC"].iloc[0]
    assert myc_row["deg_status"] == "K562_up"
    assert myc_row["regulatory_class"] == "concordant_K562"

    ctcf_a = {"chr8": pd.Series([130_300, 150_200]).to_numpy()}
    ctcf_b = {"chr8": pd.Series([425_300]).to_numpy()}
    dist_check = nearest_peak_distance(pd.Series(["chr8"]), pd.Series([130_300]), ctcf_a)
    assert dist_check[0] == 0.0

    annotated, stats_df, candidates = annotate_and_test_ctcf(
        links_merged, gene_summary, ctcf_a, ctcf_b, "K562", "GM12878", ctcf_window_bp=1000
    )
    assert not stats_df.empty
    assert "fdr_bh" in stats_df.columns
    assert candidates.iloc[0]["gene_name"] in {"MYC", "BLK"}

    # Verify all 4 plotting functions render non-empty PNGs
    plot_atac_rna_correlation(gene_summary, "K562", "GM12878", ["MYC", "BLK"], tmp_path / "corr.png")
    plot_ctcf_enrichment(stats_df, tmp_path / "ctcf.png")
    plot_distance_to_ctcf(annotated, "K562", "GM12878", tmp_path / "dist.png")
    plot_top_candidates(candidates, "K562", "GM12878", tmp_path / "cand.png")
    for fname in ("corr.png", "ctcf.png", "dist.png", "cand.png"):
        assert (tmp_path / fname).stat().st_size > 1000

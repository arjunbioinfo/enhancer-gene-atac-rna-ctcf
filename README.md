# enhancer-gene-atac-rna-ctcf
[[https://github.com/arjunbioinfo/enhancer-gene-atac-rna-ctcf/actions/workflows/ci.yml/badge.svg]](https://github.com/arjunbioinfo/enhancer-gene-atac-rna-ctcf/actions/workflows/ci.yml)
[[https://img.shields.io/badge/License-MIT-blue.svg]](LICENSE)

A reproducible, config-driven **Snakemake + Python** workflow for integrating **ATAC-seq** chromatin accessibility, **RNA-seq** gene expression, and **CTCF ChIP-seq** occupancy to link distal enhancers to target genes, test CTCF regulatory enrichment, and prioritize candidate loci.

Every push runs unit tests and executes the full Snakemake workflow on real public **ENCODE GRCh38** data (**K562** vs **GM12878** on `chr8`, including the **8q24.21 MYC / PVT1** regulatory locus) in GitHub Actions CI.

## Pipeline Overview

1. **Download (`download_encode.py`)**: Queries the ENCODE REST API for released GRCh38 ATAC-seq narrowPeak, CTCF ChIP-seq narrowPeak, and polyA+ RNA-seq gene quantification files, downloads GENCODE v29 basic GTF, and records exact `ENCFF` accessions in `data/raw/download_manifest.json`.
2. **TSS Annotation (`parse_gencode_tss.py`)**: Extracts strand-aware gene transcription start sites (protein-coding, lincRNA, antisense) from GENCODE v29.
3. **ATAC + RNA Integration (`integrate_atac_rna.py`)**: Merges condition-specific ATAC peaks into a unified consensus set, computes differential accessibility, classifies promoters vs. distal enhancers (up to 250 kb), links peaks to target genes with exponential distance decay, and classifies per-gene concordance (`concordant_K562`, `concordant_GM12878`, `chromatin_primed`, `discordant`, `baseline`).
4. **CTCF Enrichment & Prioritization (`ctcf_enrichment.py`)**: Measures exact distance to nearest CTCF ChIP-seq peak, runs two-sided Fisher's exact tests (with 95% CI and Benjamini-Hochberg FDR) and Mann-Whitney distance tests, provides a JASPAR `MA0139.1` CTCF log-odds PWM scorer, and ranks multi-omic candidate genes.
5. **Visualization (`plot_summary.py`)**: Produces four 300 DPI summary figures in `figures/`.

## Public Benchmark Data (ENCODE GRCh38)

| Cell Line | ATAC-seq | CTCF ChIP-seq | RNA-seq (polyA+) |
| :--- | :--- | :--- | :--- |
| **K562** (myeloid leukemia) | [ENCSR483RKN](https://www.encodeproject.org/experiments/ENCSR483RKN/) | [ENCSR000DWE](https://www.encodeproject.org/experiments/ENCSR000DWE/) | [ENCSR000AEM](https://www.encodeproject.org/experiments/ENCSR000AEM/) |
| **GM12878** (lymphoblastoid) | [ENCSR095QNB](https://www.encodeproject.org/experiments/ENCSR095QNB/) | [ENCSR000DZN](https://www.encodeproject.org/experiments/ENCSR000DZN/) | [ENCSR000AED](https://www.encodeproject.org/experiments/ENCSR000AED/) |

Gene TSS annotations are parsed from GENCODE release 29 basic annotation (GRCh38) to match ENCODE V29 gene quantifications.

## Example Outputs (ENCODE GRCh38 `chr8`: K562 vs. GM12878)

| Distal Enhancer Accessibility vs. Expression | CTCF Binding Enrichment |
| :---: | :---: |
| [figures/atac_rna_correlation.png] | [figures/ctcf_enrichment.png] |

| Distance from Distal Enhancers to CTCF | Top Prioritized Candidates |
| :---: | :---: |
| [figures/distance_to_ctcf.png] | [figures/top_candidates_summary.png] |

## Quick Start

```bash
# 1. Create environment (Conda)
conda env create -f environment.yml
conda activate enhancer-atac-rna-ctcf

# ...or pip:
pip install -r requirements.txt

# 2. Run unit tests (offline)
make test

# 3. Fast end-to-end run on real ENCODE chr8 data (K562 vs GM12878)
make ci

# 4. Full genome-wide run (chr1-chr22, chrX)
make run

# Or run via Docker (mounts local results/ and figures/)
docker build -t enhancer-atac-rna-ctcf .
docker run --rm -v "$(pwd)/results:/app/results" -v "$(pwd)/figures:/app/figures" enhancer-atac-rna-ctcf
```

## Outputs

* `results/enhancer_gene_links.tsv` — Consensus ATAC peaks linked to GENCODE TSSs (promoter vs. distal enhancer) with distance-decayed regulatory weights and RNA-seq log2FC.
* `results/gene_regulatory_summary.tsv` — Per-gene regulatory classification (`concordant_K562`, `concordant_GM12878`, `chromatin_primed`, `discordant`, `baseline`).
* `results/ctcf_enrichment_stats.tsv` — Two-sided Fisher's exact tests (odds ratios, 95% CI, BH-FDR) and Mann-Whitney distance tests for CTCF enrichment at condition-specific enhancers and DEG-linked enhancers.
* `results/prioritized_candidates.tsv` — Ranked multi-omic candidate genes combining concordant expression change, distal enhancer accessibility, and proximal CTCF occupancy.
* `figures/` — Publication-ready 300 DPI PNGs: `atac_rna_correlation.png`, `ctcf_enrichment.png`, `distance_to_ctcf.png`, and `top_candidates_summary.png`.

## Using Your Own Data

1. Replace the `encode` experiment accessions in `config/config.yaml` with any pair of released ENCODE GRCh38 experiments, or place local files directly in `data/raw/` named `<COND>_atac.narrowPeak.gz`, `<COND>_ctcf.narrowPeak.gz`, and `<COND>_rna_quant.tsv` (columns: `gene_id`, `TPM`).
2. Adjust genomic windows and thresholds in `config/config.yaml` (`promoter_window_bp`, `enhancer_max_distance_bp`, `ctcf_proximity_bp`, `min_tpm`, `deg_log2fc_threshold`, `atac_log2fc_threshold`).
3. For sequence-level CTCF motif scoring against a custom FASTA, import `score_ctcf_pwm` (JASPAR `MA0139.1` log-odds PWM) from `workflow/scripts/ctcf_enrichment.py`.

## Repository Layout

* `Snakefile` — End-to-end Snakemake workflow
* `Makefile` — `setup`, `test`, `ci`, `run`, and `clean` shortcuts
* `config/config.yaml` — Genome-wide K562 vs GM12878 config
* `config/ci.yaml` — Fast `chr8` (MYC locus) CI config
* `workflow/scripts/` — Modular Python scripts (`download_encode.py`, `parse_gencode_tss.py`, `integrate_atac_rna.py`, `ctcf_enrichment.py`, `plot_summary.py`)
* `notebooks/walkthrough.ipynb` — Interactive exploration notebook
* `tests/test_pipeline.py` — Unit and synthetic end-to-end smoke tests
* `environment.yml` / `requirements.txt` — Pinned Conda and pip environments
* `Dockerfile` / `.dockerignore` — Container image definition for Docker and Apptainer
* `CITATION.cff` — Software citation metadata for GitHub and Zenodo archiving
* `.github/workflows/ci.yml` — GitHub Actions CI (pytest + `chr8` ENCODE run)

## License

MIT — see [LICENSE](LICENSE).

## Author

Mallikarjuna Thippana, PhD — Computational Biologist ([arjunbioinfo.github.io](https://arjunbioinfo.github.io/)).

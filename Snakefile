"""Reproducible Snakemake workflow for ATAC-seq + RNA-seq + CTCF integration."""

configfile: "config/config.yaml"

RESULTS_DIR = config.get("output_dirs", {}).get("results", "results")
FIGURES_DIR = config.get("output_dirs", {}).get("figures", "figures")
DATA_DIR = config.get("output_dirs", {}).get("data", "data/raw")

COND_A = config["comparison"]["condition_a"]
COND_B = config["comparison"]["condition_b"]


rule all:
    input:
        f"{RESULTS_DIR}/enhancer_gene_links.tsv",
        f"{RESULTS_DIR}/gene_regulatory_summary.tsv",
        f"{RESULTS_DIR}/ctcf_enrichment_stats.tsv",
        f"{RESULTS_DIR}/prioritized_candidates.tsv",
        f"{FIGURES_DIR}/atac_rna_correlation.png",
        f"{FIGURES_DIR}/ctcf_enrichment.png",
        f"{FIGURES_DIR}/distance_to_ctcf.png",
        f"{FIGURES_DIR}/top_candidates_summary.png",


rule download_encode_data:
    output:
        atac_a=f"{DATA_DIR}/{COND_A}_atac.narrowPeak.gz",
        atac_b=f"{DATA_DIR}/{COND_B}_atac.narrowPeak.gz",
        ctcf_a=f"{DATA_DIR}/{COND_A}_ctcf.narrowPeak.gz",
        ctcf_b=f"{DATA_DIR}/{COND_B}_ctcf.narrowPeak.gz",
        rna_a=f"{DATA_DIR}/{COND_A}_rna_quant.tsv",
        rna_b=f"{DATA_DIR}/{COND_B}_rna_quant.tsv",
        gtf=f"{DATA_DIR}/gencode.basic.annotation.gtf.gz",
        manifest=f"{DATA_DIR}/download_manifest.json",
    params:
        config_path=config.get("_config_path", "config/config.yaml"),
    shell:
        """
        python workflow/scripts/download_encode.py \
            --config {params.config_path} \
            --out-dir {DATA_DIR}
        """


rule parse_gencode_tss:
    input:
        gtf=rules.download_encode_data.output.gtf,
    output:
        tss=f"{RESULTS_DIR}/gencode_tss.tsv",
    params:
        chroms=",".join(config.get("chromosomes", [])),
        biotypes=",".join(config.get("annotation", {}).get("biotypes", ["protein_coding", "lincRNA", "antisense"])),
    shell:
        """
        python workflow/scripts/parse_gencode_tss.py \
            --gtf {input.gtf} \
            --out {output.tss} \
            --chromosomes "{params.chroms}" \
            --biotypes "{params.biotypes}"
        """


rule integrate_atac_rna:
    input:
        atac_a=rules.download_encode_data.output.atac_a,
        atac_b=rules.download_encode_data.output.atac_b,
        rna_a=rules.download_encode_data.output.rna_a,
        rna_b=rules.download_encode_data.output.rna_b,
        tss=rules.parse_gencode_tss.output.tss,
    output:
        links=f"{RESULTS_DIR}/enhancer_gene_links.tsv",
        gene_summary=f"{RESULTS_DIR}/gene_regulatory_summary.tsv",
    params:
        cond_a=COND_A,
        cond_b=COND_B,
        chroms=",".join(config.get("chromosomes", [])),
        promoter_bp=config["parameters"]["promoter_window_bp"],
        enhancer_max_bp=config["parameters"]["enhancer_max_distance_bp"],
        min_tpm=config["parameters"]["min_tpm"],
        deg_lfc=config["parameters"]["deg_log2fc_threshold"],
        atac_lfc=config["parameters"]["atac_log2fc_threshold"],
    shell:
        """
        python workflow/scripts/integrate_atac_rna.py \
            --atac-a {input.atac_a} \
            --atac-b {input.atac_b} \
            --rna-a {input.rna_a} \
            --rna-b {input.rna_b} \
            --tss {input.tss} \
            --cond-a {params.cond_a} \
            --cond-b {params.cond_b} \
            --chromosomes "{params.chroms}" \
            --promoter-bp {params.promoter_bp} \
            --enhancer-max-bp {params.enhancer_max_bp} \
            --min-tpm {params.min_tpm} \
            --deg-lfc {params.deg_lfc} \
            --atac-lfc {params.atac_lfc} \
            --out-links {output.links} \
            --out-genes {output.gene_summary}
        """


rule ctcf_enrichment:
    input:
        links=rules.integrate_atac_rna.output.links,
        gene_summary=rules.integrate_atac_rna.output.gene_summary,
        ctcf_a=rules.download_encode_data.output.ctcf_a,
        ctcf_b=rules.download_encode_data.output.ctcf_b,
    output:
        stats=f"{RESULTS_DIR}/ctcf_enrichment_stats.tsv",
        annotated_links=f"{RESULTS_DIR}/enhancer_gene_ctcf_annotated.tsv",
        candidates=f"{RESULTS_DIR}/prioritized_candidates.tsv",
    params:
        cond_a=COND_A,
        cond_b=COND_B,
        chroms=",".join(config.get("chromosomes", [])),
        ctcf_window_bp=config["parameters"]["ctcf_proximity_bp"],
        seed=config["parameters"]["random_seed"],
    shell:
        """
        python workflow/scripts/ctcf_enrichment.py \
            --links {input.links} \
            --genes {input.gene_summary} \
            --ctcf-a {input.ctcf_a} \
            --ctcf-b {input.ctcf_b} \
            --cond-a {params.cond_a} \
            --cond-b {params.cond_b} \
            --chromosomes "{params.chroms}" \
            --ctcf-window-bp {params.ctcf_window_bp} \
            --seed {params.seed} \
            --out-stats {output.stats} \
            --out-links {output.annotated_links} \
            --out-candidates {output.candidates}
        """


rule plot_summary:
    input:
        links=rules.ctcf_enrichment.output.annotated_links,
        genes=rules.integrate_atac_rna.output.gene_summary,
        stats=rules.ctcf_enrichment.output.stats,
        candidates=rules.ctcf_enrichment.output.candidates,
    output:
        corr_fig=f"{FIGURES_DIR}/atac_rna_correlation.png",
        ctcf_fig=f"{FIGURES_DIR}/ctcf_enrichment.png",
        dist_fig=f"{FIGURES_DIR}/distance_to_ctcf.png",
        cand_fig=f"{FIGURES_DIR}/top_candidates_summary.png",
    params:
        cond_a=COND_A,
        cond_b=COND_B,
        highlight_genes=",".join(config.get("plotting", {}).get("highlight_genes", ["MYC", "GATA1", "PVT1", "CD79A", "PAX5"])),
    shell:
        """
        python workflow/scripts/plot_summary.py \
            --links {input.links} \
            --genes {input.genes} \
            --stats {input.stats} \
            --candidates {input.candidates} \
            --cond-a {params.cond_a} \
            --cond-b {params.cond_b} \
            --highlight-genes "{params.highlight_genes}" \
            --out-dir {FIGURES_DIR}
        """

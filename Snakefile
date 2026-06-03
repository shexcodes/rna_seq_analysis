# =============================================================================
#  RNA-Seq analysis pipeline
#  Paired-end reads -> QC -> trimming -> alignment -> counts -> DE analysis
#
#  Steps:
#    1. FastQC            quality control of raw reads
#    2. flexbar           adapter & quality trimming
#    2b. FastQC           quality control of trimmed reads
#    3. STAR              build genome index (once) + align reads
#    4. featureCounts     count reads per gene
#    5. DESeq2 (R)        differential expression + plots
#    6. MultiQC           aggregate all QC into one report
#
#  Run from the repository root:
#      snakemake -n                 # dry run (plan only)
#      snakemake --cores 4          # run everything
# =============================================================================

import pandas as pd

configfile: "config/config.yaml"

# ---- Sample sheet -----------------------------------------------------------
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
    .set_index("sample", drop=False)
    .sort_index()
)
SAMPLES = samples["sample"].tolist()
READS = ["1", "2"]


wildcard_constraints:
    sample="|".join(SAMPLES),
    read="|".join(READS),


def fq(sample, read):
    """Return the raw FASTQ path for a given sample and read (1 or 2)."""
    return samples.loc[sample, "fq1" if read == "1" else "fq2"]


# ---- Target rule: everything the pipeline should produce --------------------
rule all:
    input:
        "results/deseq2/results.tsv",
        "results/deseq2/normalized_counts.tsv",
        "results/deseq2/plots/pca.png",
        "results/deseq2/plots/ma_plot.png",
        "results/deseq2/plots/volcano.png",
        "results/deseq2/plots/sample_distances.png",
        "results/deseq2/plots/top_genes_heatmap.png",
        "results/enrichment/go_enrichment.tsv",
        "results/enrichment/kegg_enrichment.tsv",
        "results/enrichment/plots/go_dotplot.png",
        "results/enrichment/plots/kegg_dotplot.png",
        "results/qc/multiqc_report.html",


# ---- 1. Quality control of raw reads ----------------------------------------
rule fastqc_raw:
    input:
        lambda wc: fq(wc.sample, wc.read),
    output:
        html="results/qc/fastqc_raw/{sample}_R{read}_fastqc.html",
        zip="results/qc/fastqc_raw/{sample}_R{read}_fastqc.zip",
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/fastqc_raw/{sample}_R{read}.log",
    threads: 1
    shell:
        # FastQC names outputs after the input file, so run in a temp dir
        # and rename to deterministic, wildcard-driven names.
        r"""
        tmp=$(mktemp -d)
        fastqc -q -t {threads} -o "$tmp" {input} > {log} 2>&1
        mv "$tmp"/*_fastqc.html {output.html}
        mv "$tmp"/*_fastqc.zip  {output.zip}
        rm -rf "$tmp"
        """


# ---- 2. Adapter & quality trimming (flexbar) --------------------------------
rule flexbar_trim:
    input:
        r1=lambda wc: fq(wc.sample, "1"),
        r2=lambda wc: fq(wc.sample, "2"),
        adapters=config["reference"]["adapters"],
    output:
        r1="results/trimmed/{sample}_1.fastq.gz",
        r2="results/trimmed/{sample}_2.fastq.gz",
    params:
        prefix="results/trimmed/{sample}",
        pre_trim_left=config["params"]["flexbar"]["pre_trim_left"],
        min_overlap=config["params"]["flexbar"]["adapter_min_overlap"],
        min_len=config["params"]["flexbar"]["min_read_length"],
    conda:
        "envs/flexbar.yaml"
    log:
        "logs/flexbar/{sample}.log",
    threads: config["threads"]["flexbar"]
    shell:
        # --pre-trim-left (-x) removes the first N bases of every read.
        # --zip-output GZ writes gzipped FASTQ so STAR can read it with zcat.
        r"""
        flexbar \
            --reads {input.r1} --reads2 {input.r2} \
            --adapters {input.adapters} --adapters2 {input.adapters} \
            --pre-trim-left {params.pre_trim_left} \
            --adapter-min-overlap {params.min_overlap} \
            --min-read-length {params.min_len} \
            --zip-output GZ \
            --threads {threads} \
            --target {params.prefix} > {log} 2>&1
        """


# ---- 2b. Quality control of trimmed reads -----------------------------------
rule fastqc_trimmed:
    input:
        "results/trimmed/{sample}_{read}.fastq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}_R{read}_fastqc.html",
        zip="results/qc/fastqc_trimmed/{sample}_R{read}_fastqc.zip",
    conda:
        "envs/fastqc.yaml"
    log:
        "logs/fastqc_trimmed/{sample}_R{read}.log",
    threads: 1
    shell:
        r"""
        tmp=$(mktemp -d)
        fastqc -q -t {threads} -o "$tmp" {input} > {log} 2>&1
        mv "$tmp"/*_fastqc.html {output.html}
        mv "$tmp"/*_fastqc.zip  {output.zip}
        rm -rf "$tmp"
        """


# ---- 3a. STAR genome index (built once, reused by every sample) -------------
rule star_index:
    input:
        genome=config["reference"]["genome"],
        gtf=config["reference"]["annotation"],
    output:
        directory("results/star_index"),
    params:
        sa_index_nbases=config["params"]["star"]["sa_index_nbases"],
        sjdb_overhang=config["params"]["star"]["sjdb_overhang"],
    conda:
        "envs/star.yaml"
    log:
        "logs/star/index.log",
    threads: config["threads"]["star"]
    shell:
        r"""
        mkdir -p {output}
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            --genomeSAindexNbases {params.sa_index_nbases} \
            --runThreadN {threads} > {log} 2>&1
        """


# ---- 3b. STAR alignment (per sample, outputs a sorted BAM) ------------------
rule star_align:
    input:
        r1="results/trimmed/{sample}_1.fastq.gz",
        r2="results/trimmed/{sample}_2.fastq.gz",
        index="results/star_index",
    output:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log_final="results/star/{sample}/Log.final.out",
    params:
        prefix="results/star/{sample}/",
    conda:
        "envs/star.yaml"
    log:
        "logs/star/{sample}.log",
    threads: config["threads"]["star"]
    shell:
        r"""
        STAR --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} > {log} 2>&1
        """


# ---- 4. Feature counting (one featureCounts call over all BAMs) -------------
rule feature_counts:
    input:
        bams=expand(
            "results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES
        ),
        gtf=config["reference"]["annotation"],
    output:
        counts="results/counts/counts.tsv",
        summary="results/counts/counts.tsv.summary",
    params:
        feature_type=config["params"]["featurecounts"]["feature_type"],
        attribute=config["params"]["featurecounts"]["attribute"],
    conda:
        "envs/subread.yaml"
    log:
        "logs/featurecounts.log",
    threads: config["threads"]["featurecounts"]
    shell:
        # -p with --countReadPairs counts fragments (read pairs) for PE data.
        r"""
        featureCounts \
            -p --countReadPairs \
            -T {threads} \
            -t {params.feature_type} \
            -g {params.attribute} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} > {log} 2>&1
        """


# ---- 5. Differential expression analysis (DESeq2) ---------------------------
rule deseq2:
    input:
        counts="results/counts/counts.tsv",
        samples=config["samples"],
    output:
        results="results/deseq2/results.tsv",
        normalized="results/deseq2/normalized_counts.tsv",
        pca="results/deseq2/plots/pca.png",
        ma="results/deseq2/plots/ma_plot.png",
        volcano="results/deseq2/plots/volcano.png",
        distances="results/deseq2/plots/sample_distances.png",
        heatmap="results/deseq2/plots/top_genes_heatmap.png",
    params:
        contrast=config["deseq2"]["contrast"],
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
        top_n=config["deseq2"]["top_n_heatmap"],
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2.log",
    script:
        "scripts/deseq2.R"


# ---- 5b. Functional enrichment of significant genes (GO / KEGG) -------------
rule enrichment:
    input:
        results="results/deseq2/results.tsv",
    output:
        go_table="results/enrichment/go_enrichment.tsv",
        kegg_table="results/enrichment/kegg_enrichment.tsv",
        go_plot="results/enrichment/plots/go_dotplot.png",
        kegg_plot="results/enrichment/plots/kegg_dotplot.png",
    params:
        orgdb=config["enrichment"]["orgdb"],
        kegg_organism=config["enrichment"]["kegg_organism"],
        gene_id_type=config["enrichment"]["gene_id_type"],
        go_ontology=config["enrichment"]["go_ontology"],
        qvalue_cutoff=config["enrichment"]["qvalue_cutoff"],
        show_categories=config["enrichment"]["show_categories"],
        # significance thresholds reused from the DE step for consistency
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
    conda:
        "envs/enrichment.yaml"
    log:
        "logs/enrichment.log",
    script:
        "scripts/enrichment.R"


# ---- 6. Aggregate QC report (MultiQC) ---------------------------------------
rule multiqc:
    input:
        expand(
            "results/qc/fastqc_raw/{sample}_R{read}_fastqc.zip",
            sample=SAMPLES,
            read=READS,
        ),
        expand(
            "results/qc/fastqc_trimmed/{sample}_R{read}_fastqc.zip",
            sample=SAMPLES,
            read=READS,
        ),
        expand("logs/flexbar/{sample}.log", sample=SAMPLES),
        expand("results/star/{sample}/Log.final.out", sample=SAMPLES),
        "results/counts/counts.tsv.summary",
    output:
        "results/qc/multiqc_report.html",
    conda:
        "envs/multiqc.yaml"
    log:
        "logs/multiqc.log",
    shell:
        r"""
        multiqc results/ logs/ \
            --force \
            --outdir results/qc \
            --filename multiqc_report.html > {log} 2>&1
        """

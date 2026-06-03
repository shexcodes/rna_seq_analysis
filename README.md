# RNA-Seq Analysis Notebook

Downstream RNA-Seq analysis in Python — from a gene-level **count matrix** to
**differential expression** and **GO/KEGG enrichment**, with every step shown
and executed in one notebook.

This is the interactive companion to the `rnaseq-snakemake` pipeline, which
produces the count matrix from raw reads (FastQC → flexbar → STAR →
featureCounts). The two together cover a complete RNA-Seq workflow.

## Contents

```
rnaseq-analysis-notebook/
├── notebooks/
│   └── rnaseq_analysis.ipynb   # the analysis (pre-run, outputs embedded)
├── data/
│   ├── counts.tsv              # gene x sample count matrix
│   ├── sample_metadata.csv     # sample, condition, replicate
│   └── generate_counts.py      # reproducible data generator (seed=42)
├── results/                    # written by the notebook (git-ignored)
├── requirements.txt
└── .gitignore
```

## What the notebook does

1. Documents the upstream (reads → counts) commands for reference
2. Loads the count matrix and sample metadata
3. Exploratory QC (library sizes, detected genes, count distributions)
4. Filters low-count genes
5. Differential expression with **PyDESeq2** (G2 vs G1)
6. Normalization (logCPM), **PCA**, sample-correlation heatmap
7. **MA** and **volcano** plots
8. Heatmap of the top differentially expressed genes
9. **GO / KEGG** over-representation with **gseapy** (Enrichr)
10. Summary and saved outputs

## Quickstart

```bash
pip install -r requirements.txt

# (optional) regenerate the synthetic dataset
python data/generate_counts.py

jupyter lab notebooks/rnaseq_analysis.ipynb
```

## Using your own data

Replace `data/counts.tsv` (genes × samples, tab-separated, first column `gene`)
and `data/sample_metadata.csv` (`sample`, `condition`, ...). If your conditions
aren't named `G1`/`G2`, update the contrast and colors in the setup cell.

> The included dataset is **simulated** (two conditions × three replicates, with
> real human cell-cycle / interferon / apoptosis genes made differentially
> expressed) so the notebook runs without large reference files. The GO/KEGG
> step queries Enrichr and therefore needs internet access.

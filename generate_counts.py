"""
Generate a synthetic but realistic bulk RNA-seq count matrix + sample metadata.

Design: two conditions (G1 = control, G2 = treatment), three replicates each.
A curated set of *real* human genes from coherent pathways (cell cycle, type-I
interferon / inflammation, apoptosis) are made differentially expressed so that
downstream differential-expression and GO/KEGG enrichment recover sensible
biology. All other genes are background with no condition effect. Counts are
drawn from a negative-binomial model with per-gene dispersion and per-sample
size factors, mimicking real RNA-seq count behaviour.

Reproducible (fixed seed). Outputs (written next to this script):
    counts.tsv            gene x sample integer count matrix (tab-separated)
    sample_metadata.csv   sample, condition, replicate
"""

from pathlib import Path
import numpy as np
import pandas as pd

RNG = np.random.default_rng(42)
OUT = Path(__file__).resolve().parent

# ---- Samples ---------------------------------------------------------------
samples = ["G1_rep1", "G1_rep2", "G1_rep3", "G2_rep1", "G2_rep2", "G2_rep3"]
condition = ["G1", "G1", "G1", "G2", "G2", "G2"]
replicate = [1, 2, 3, 1, 2, 3]

# ---- Real human genes with a coherent expected signal ----------------------
# Up in G2 (proliferation + interferon/inflammation)
UP_GENES = [
    # cell cycle / mitosis
    "CDK1", "CCNB1", "CCNB2", "CCNA2", "CCNE1", "CDC20", "CDC25A", "CDC25C",
    "BUB1", "BUB1B", "MAD2L1", "PLK1", "AURKA", "AURKB", "TOP2A", "MKI67",
    "FOXM1", "KIF11", "KIF23", "CENPA", "CENPE", "CENPF", "NUSAP1", "TPX2",
    "UBE2C", "BIRC5", "TTK", "ESPL1", "PTTG1", "RRM2", "TYMS", "PCNA",
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "E2F1",
    # type-I interferon / inflammation
    "IL6", "TNF", "CXCL8", "CXCL10", "CCL2", "CCL5", "STAT1", "STAT2",
    "IRF1", "IRF7", "NFKB1", "IL1B", "ISG15", "OAS1", "OAS2", "MX1",
    "IFIT1", "IFIT3", "GBP1",
]
# Down in G2 (apoptosis)
DOWN_GENES = [
    "BAX", "BAK1", "CASP3", "CASP7", "CASP8", "CASP9", "TP53", "BBC3",
    "PMAIP1", "APAF1", "CYCS", "FAS", "FADD", "BID", "DIABLO",
]

N_GENES = 2000
n_bg = N_GENES - len(UP_GENES) - len(DOWN_GENES)
bg_genes = [f"BG{i:05d}" for i in range(1, n_bg + 1)]
genes = UP_GENES + DOWN_GENES + bg_genes
assert len(genes) == len(set(genes)) == N_GENES

# ---- Per-gene baseline expression and dispersion ---------------------------
base_mean = np.exp(RNG.normal(3.5, 2.0, N_GENES)).clip(5, 50000)
disp = RNG.uniform(0.10, 0.40, N_GENES)              # NB dispersion
size_factors = RNG.uniform(0.8, 1.2, len(samples))   # per-sample depth

# ---- True log2 fold changes (G2 vs G1) -------------------------------------
lfc = np.zeros(N_GENES)
lfc[: len(UP_GENES)] = RNG.uniform(1.5, 3.0, len(UP_GENES))
lfc[len(UP_GENES) : len(UP_GENES) + len(DOWN_GENES)] = -RNG.uniform(
    1.5, 3.0, len(DOWN_GENES)
)

# ---- Simulate counts (negative binomial) -----------------------------------
counts = np.zeros((N_GENES, len(samples)), dtype=int)
nb_n = 1.0 / disp
for j, (cond, sf) in enumerate(zip(condition, size_factors)):
    fold = 2.0 ** lfc if cond == "G2" else np.ones(N_GENES)
    mean_j = base_mean * fold * sf
    p = nb_n / (nb_n + mean_j)
    counts[:, j] = RNG.negative_binomial(nb_n, p)

# ---- Write outputs ----------------------------------------------------------
counts_df = pd.DataFrame(counts, index=genes, columns=samples)
counts_df.index.name = "gene"
counts_df.to_csv(OUT / "counts.tsv", sep="\t")

meta = pd.DataFrame(
    {"sample": samples, "condition": condition, "replicate": replicate}
)
meta.to_csv(OUT / "sample_metadata.csv", index=False)

print(f"Wrote {OUT/'counts.tsv'}  ({counts_df.shape[0]} genes x {counts_df.shape[1]} samples)")
print(f"Wrote {OUT/'sample_metadata.csv'}")
print(f"Differentially expressed (true): {len(UP_GENES)} up, {len(DOWN_GENES)} down")

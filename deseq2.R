# =============================================================================
#  Differential expression analysis with DESeq2
#  Run by Snakemake via the `script:` directive (inputs/outputs/params are
#  available in the `snakemake` S4 object). It reads the featureCounts matrix,
#  fits a DESeq2 model, and writes a results table, normalized counts, and a
#  standard set of diagnostic plots.
# =============================================================================

# Redirect all R messages/output to the Snakemake log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# ---- Inputs & parameters ----------------------------------------------------
counts_file  <- snakemake@input[["counts"]]
samples_file <- snakemake@input[["samples"]]
contrast     <- unlist(snakemake@params[["contrast"]])   # c(factor, test, ref)
alpha        <- as.numeric(snakemake@params[["alpha"]])
lfc_thresh   <- as.numeric(snakemake@params[["lfc_threshold"]])
top_n        <- as.integer(snakemake@params[["top_n"]])

factor_name <- contrast[1]
level_test  <- contrast[2]
level_ref   <- contrast[3]

# ---- Read the featureCounts matrix ------------------------------------------
# featureCounts output: a comment line (#), then a header with columns
# Geneid, Chr, Start, End, Strand, Length, then one column per BAM.
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
count_cols <- 7:ncol(fc)
mat <- as.matrix(fc[, count_cols])
rownames(mat) <- fc[["Geneid"]]

# featureCounts labels count columns with the BAM path
# (results/star/<sample>/Aligned.sortedByCoord.out.bam); recover the sample id.
colnames(mat) <- basename(dirname(colnames(fc)[count_cols]))
storage.mode(mat) <- "integer"

# ---- Sample metadata --------------------------------------------------------
coldata <- read.delim(samples_file, comment.char = "#", check.names = FALSE,
                      stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(mat), , drop = FALSE]    # align to matrix order
coldata[[factor_name]] <- relevel(factor(coldata[[factor_name]]),
                                  ref = level_ref)

stopifnot(all(rownames(coldata) == colnames(mat)))
message("Loaded ", nrow(mat), " genes across ", ncol(mat), " samples.")

# ---- Build DESeq2 dataset & run ---------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData   = coldata,
  design    = as.formula(paste("~", factor_name))
)

# Pre-filter: keep genes with >=10 counts in at least the size of the
# smallest group (removes noise, speeds up and stabilises fitting).
min_group <- min(table(coldata[[factor_name]]))
keep <- rowSums(counts(dds) >= 10) >= min_group
dds  <- dds[keep, ]
message(sum(keep), " genes pass the low-count filter.")

dds <- DESeq(dds)

# ---- Results (with log2 fold-change shrinkage) ------------------------------
res <- results(dds, contrast = c(factor_name, level_test, level_ref),
               alpha = alpha)

coef_name <- paste0(factor_name, "_", level_test, "_vs_", level_ref)
res_shrunk <- tryCatch(
  lfcShrink(dds, coef = coef_name, type = "apeglm"),
  error = function(e) {
    message("apeglm shrinkage unavailable (", conditionMessage(e),
            "); falling back to normal.")
    lfcShrink(dds, contrast = c(factor_name, level_test, level_ref),
              type = "normal")
  }
)

res_df <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), c("gene", setdiff(names(res_df), "gene"))]
write.table(res_df, snakemake@output[["results"]],
            sep = "\t", quote = FALSE, row.names = FALSE)

n_sig <- sum(res_df$padj < alpha & abs(res_df$log2FoldChange) > lfc_thresh,
             na.rm = TRUE)
message(n_sig, " genes significant at padj < ", alpha,
        " and |log2FC| > ", lfc_thresh, ".")

# ---- Normalized counts ------------------------------------------------------
norm <- counts(dds, normalized = TRUE)
write.table(
  data.frame(gene = rownames(norm), norm, check.names = FALSE),
  snakemake@output[["normalized"]],
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---- Variance-stabilising transform for visualisation -----------------------
vsd <- tryCatch(
  vst(dds, blind = FALSE),
  error = function(e) varianceStabilizingTransformation(dds, blind = FALSE)
)

# ---- PCA --------------------------------------------------------------------
pca <- plotPCA(vsd, intgroup = factor_name, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, color = .data[[factor_name]], label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", pv[1], "% variance")) +
  ylab(paste0("PC2: ", pv[2], "% variance")) +
  labs(title = "PCA of samples (VST)", color = factor_name) +
  theme_bw()
ggsave(snakemake@output[["pca"]], p_pca, width = 7, height = 5, dpi = 150)

# ---- MA plot ----------------------------------------------------------------
png(snakemake@output[["ma"]], width = 1000, height = 800, res = 150)
plotMA(res_shrunk, alpha = alpha,
       main = paste("MA plot:", level_test, "vs", level_ref))
dev.off()

# ---- Volcano plot -----------------------------------------------------------
v <- res_df
v$padj[is.na(v$padj)] <- 1
v$sig <- with(v, ifelse(
  padj < alpha & abs(log2FoldChange) > lfc_thresh,
  ifelse(log2FoldChange > 0, "Up", "Down"),
  "n.s."
))
p_volcano <- ggplot(v, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c(Up = "#d73027", Down = "#4575b4",
                                "n.s." = "grey70")) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  labs(title = paste("Volcano:", level_test, "vs", level_ref),
       x = "log2 fold-change", y = "-log10 adjusted p", color = NULL) +
  theme_bw()
ggsave(snakemake@output[["volcano"]], p_volcano, width = 7, height = 5, dpi = 150)

# ---- Sample-to-sample distance heatmap --------------------------------------
sampleDists <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
png(snakemake@output[["distances"]], width = 900, height = 800, res = 150)
pheatmap(dist_mat,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-sample distances (VST)")
dev.off()

# ---- Heatmap of the top differential genes ----------------------------------
sel <- head(order(res$padj), top_n)
ann <- data.frame(coldata[[factor_name]], row.names = rownames(coldata))
colnames(ann) <- factor_name
png(snakemake@output[["heatmap"]], width = 900, height = 1000, res = 150)
pheatmap(assay(vsd)[sel, , drop = FALSE],
         scale = "row",
         annotation_col = ann,
         show_rownames = TRUE,
         main = paste("Top", length(sel), "genes by adjusted p"))
dev.off()

message("DESeq2 analysis complete.")

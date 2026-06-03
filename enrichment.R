# =============================================================================
#  GO & KEGG over-representation analysis (ORA) on significant DE genes
#  Run by Snakemake via the `script:` directive.
#
#  Takes the DESeq2 results table, selects significant genes (same thresholds
#  as the DE step), and tests whether GO terms / KEGG pathways are enriched
#  among them relative to all tested genes (the "universe"). Always writes
#  every declared output, even when no terms pass — so the pipeline never
#  fails on a sparse dataset.
# =============================================================================

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

suppressMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

# ---- Inputs & parameters ----------------------------------------------------
res_file   <- snakemake@input[["results"]]
orgdb_name <- snakemake@params[["orgdb"]]
kegg_org   <- snakemake@params[["kegg_organism"]]
id_type    <- snakemake@params[["gene_id_type"]]   # keytype of input IDs
go_ont     <- snakemake@params[["go_ontology"]]    # BP / MF / CC / ALL
alpha      <- as.numeric(snakemake@params[["alpha"]])
lfc_thresh <- as.numeric(snakemake@params[["lfc_threshold"]])
qcut       <- as.numeric(snakemake@params[["qvalue_cutoff"]])
show_n     <- as.integer(snakemake@params[["show_categories"]])

suppressMessages(library(orgdb_name, character.only = TRUE))
orgdb <- get(orgdb_name)

# ---- Small helpers so every output always exists ----------------------------
write_empty_tsv <- function(path) {
  writeLines(
    paste("ID", "Description", "GeneRatio", "BgRatio", "pvalue",
          "p.adjust", "qvalue", "geneID", "Count", sep = "\t"),
    path
  )
}
placeholder_png <- function(path, msg) {
  png(path, width = 900, height = 700, res = 150)
  par(mar = c(0, 0, 0, 0)); plot.new(); text(0.5, 0.5, msg, cex = 1.1)
  dev.off()
}
save_enrichment <- function(obj, table_path, plot_path, title) {
  df <- if (is.null(obj)) data.frame() else as.data.frame(obj)
  if (nrow(df) > 0) {
    write.table(df, table_path, sep = "\t", quote = FALSE, row.names = FALSE)
    p <- tryCatch({
      d <- dotplot(obj, showCategory = show_n) + ggtitle(title)
      if (go_ont == "ALL" && "ONTOLOGY" %in% colnames(df)) {
        d <- d + facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")
      }
      d
    }, error = function(e) {
      message("dotplot failed (", conditionMessage(e), "); using barplot.")
      barplot(obj, showCategory = show_n) + ggtitle(title)
    })
    ggsave(plot_path, p, width = 8, height = 9, dpi = 150, limitsize = FALSE)
  } else {
    message("No enriched terms for: ", title)
    write_empty_tsv(table_path)
    placeholder_png(plot_path, paste0(title, "\n(no significant terms)"))
  }
}

# ---- Read DE results & define gene sets -------------------------------------
res <- read.delim(res_file, check.names = FALSE)
strip_version <- function(x) sub("\\..*$", "", x)   # ENSG000...4 -> ENSG000...
res$gene_clean <- strip_version(res$gene)

sig <- subset(res, !is.na(padj) & padj < alpha & abs(log2FoldChange) > lfc_thresh)
sig_genes <- unique(sig$gene_clean)
universe  <- unique(res$gene_clean)
message(length(sig_genes), " significant genes; universe = ",
        length(universe), " genes.")

# ---- GO over-representation -------------------------------------------------
ego <- NULL
if (length(sig_genes) > 0) {
  ego <- tryCatch(
    enrichGO(gene = sig_genes, universe = universe, OrgDb = orgdb,
             keyType = id_type, ont = go_ont, pAdjustMethod = "BH",
             qvalueCutoff = qcut, readable = TRUE),
    error = function(e) {
      message("enrichGO failed: ", conditionMessage(e)); NULL
    }
  )
}
save_enrichment(ego,
                snakemake@output[["go_table"]],
                snakemake@output[["go_plot"]],
                "GO enrichment (significant genes)")

# ---- KEGG over-representation (needs Entrez IDs; queries KEGG online) --------
sig_entrez <- tryCatch(
  bitr(sig_genes, fromType = id_type, toType = "ENTREZID", OrgDb = orgdb),
  error = function(e) NULL
)
uni_entrez <- tryCatch(
  bitr(universe, fromType = id_type, toType = "ENTREZID", OrgDb = orgdb),
  error = function(e) NULL
)

kk <- NULL
if (!is.null(sig_entrez) && nrow(sig_entrez) > 0) {
  kk <- tryCatch(
    enrichKEGG(gene = sig_entrez$ENTREZID,
               universe = if (!is.null(uni_entrez)) uni_entrez$ENTREZID else NULL,
               organism = kegg_org, pAdjustMethod = "BH",
               qvalueCutoff = qcut),
    error = function(e) {
      message("enrichKEGG failed (KEGG needs internet access): ",
              conditionMessage(e)); NULL
    }
  )
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    kk <- setReadable(kk, orgdb, keyType = "ENTREZID")  # Entrez -> symbols
  }
}
save_enrichment(kk,
                snakemake@output[["kegg_table"]],
                snakemake@output[["kegg_plot"]],
                "KEGG enrichment (significant genes)")

message("Enrichment analysis complete.")

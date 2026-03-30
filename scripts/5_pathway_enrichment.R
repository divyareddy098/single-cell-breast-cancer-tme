# Purpose: Perform GO enrichment on differential expression results

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
})

# Input
input_csv <- "results/tables/differential_expression_Cluster_0_vs_Cluster_1.csv"

# Create output directories
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# Load DE results
if (!file.exists(input_csv)) {
  stop("DE results file not found: ", input_csv,
       "\nRun scripts/04_differential_expression.R first.")
}

de_res <- read.csv(input_csv, row.names = 1)


# Select significant genes
de_res$gene <- rownames(de_res)

sig_genes <- de_res %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

if (length(sig_genes) == 0) {
  stop("No significant genes found for enrichment.")
}

# Convert gene symbols to ENTREZ IDs
gene_df <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) %>%
  distinct(ENTREZID, .keep_all = TRUE)

if (nrow(gene_df) == 0) {
  stop("No genes could be mapped from SYMBOL to ENTREZID.")
}

# GO enrichment analysis
ego <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# Save results
write.csv(as.data.frame(ego), "results/tables/go_enrichment_results.csv", row.names = FALSE)

# Plot enrichment results
png("results/figures/go_enrichment_dotplot.png", width = 1200, height = 800)
print(dotplot(ego, showCategory = 15) + ggtitle("GO Biological Process Enrichment"))
dev.off()

message("GO enrichment results saved.")

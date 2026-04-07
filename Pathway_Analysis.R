

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(EnhancedVolcano)
library(pathview)

# -------------------
# 1. Load and prepare DESeq2 results
# -------------------
res <- read.csv("DGESeq_results_cleaned.csv", stringsAsFactors = FALSE)
glimpse(res)
summary(res$padj)

# Filter for significant genes (padj < 0.05)
sig_genes <- res %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(padj)
cat("Number of significant genes (padj < 0.05):", nrow(sig_genes), "\n")

# -------------------
# 2. Prepare gene lists
# -------------------
# All genes with valid symbols (background)
all_genes <- res$symbol[!is.na(res$symbol) & res$symbol != ""]

# Ranked genes for GSEA
ranked_genes <- res %>%
  filter(!is.na(symbol) & symbol != "" & !is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(log2FoldChange, name = symbol)

# Up- and down-regulated gene lists
up_genes <- sig_genes %>% filter(log2FoldChange > 1) %>% pull(symbol) %>% unique()
down_genes <- sig_genes %>% filter(log2FoldChange < -1) %>% pull(symbol) %>% unique()
cat("Up-regulated genes:", length(up_genes), "\n")
cat("Down-regulated genes:", length(down_genes), "\n")

# -------------------
# 3. Prepare Entrez ID mapping for GSEA
# -------------------
gene_entrez <- data.frame(
  SYMBOL = names(ranked_genes),
  ENTREZID = mapIds(
    org.Hs.eg.db,
    keys = names(ranked_genes),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  ),
  stringsAsFactors = FALSE
)

# Filter ranked genes to only those with valid Entrez IDs
valid_symbols <- names(ranked_genes) %in% gene_entrez$SYMBOL
ranked_entrez <- ranked_genes[valid_symbols]
entrez_ids <- gene_entrez$ENTREZID[match(names(ranked_entrez), gene_entrez$SYMBOL)]
keep <- !is.na(entrez_ids)
ranked_entrez <- ranked_entrez[keep]
entrez_ids <- entrez_ids[keep]
names(ranked_entrez) <- entrez_ids

# Collapse duplicate Entrez IDs by max log2FC
ranked_entrez_df <- data.frame(EntrezID = names(ranked_entrez), Score = as.numeric(ranked_entrez))
ranked_entrez_unique <- ranked_entrez_df %>%
  group_by(EntrezID) %>%
  summarise(Score = max(Score), .groups = "drop") %>%
  arrange(desc(Score))
ranked_entrez_final <- ranked_entrez_unique$Score
names(ranked_entrez_final) <- ranked_entrez_unique$EntrezID

# -------------------
# 4. Over-Representation Analysis (ORA) - GO & KEGG
# -------------------
ego_bp <- enrichGO(
  gene          = up_genes,
  universe      = all_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

ego_mf <- enrichGO(
  gene          = up_genes,
  universe      = all_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

ekegg <- enrichKEGG(
  gene          = up_genes,
  universe      = all_genes,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# -------------------
# 5. GSEA KEGG
# -------------------
gsea_kegg <- gseKEGG(
  geneList     = ranked_entrez_final,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# -------------------
# 6. Visualization
# -------------------
# GO Dotplot
dotplot(ego_bp, showCategory = Inf, label_format = 70) +
  ggtitle("GO BP - Upregulated Genes") +
  theme(axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))

# GO Barplot
barplot(ego_bp, showCategory = Inf, label_format = 70) +
  ggtitle("GO BP Enrichment") +
  theme(axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))

# KEGG Dotplot
if (!is.null(ekegg)) {
  dotplot(ekegg, showCategory = Inf, label_format = 70) +
    ggtitle("KEGG Pathway Enrichment") +
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
}

# GSEA Ridgeplot
ridgeplot(gsea_kegg, showCategory = 10) +
  ggtitle("GSEA KEGG Pathways")

# Enhanced Volcano
EnhancedVolcano(
  res,
  lab = res$symbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Differential Expression',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0,
  colAlpha = 0.7
)

# -------------------
# 7. Save results
# -------------------
write.csv(as.data.frame(ego_bp), "GO_BP_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "KEGG_enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), "GSEA_KEGG.csv", row.names = FALSE)

if ("dplyr" %in% loadedNamespaces()) {
  sig_genes %>%
    dplyr::select(gene, symbol, baseMean, log2FoldChange, pvalue, padj) %>%
    write.csv("Significant_genes.csv", row.names = FALSE)
} else {
  write.csv(sig_genes[, c("gene","symbol","baseMean","log2FoldChange","pvalue","padj")],
            "Significant_genes.csv",
            row.names = FALSE)
}

message("All enrichment and DE results saved successfully.")
# ============================================================
# Script: 06_enrichment.R
# Purpose: GO and KEGG pathway enrichment analysis on ME11
#          genes to characterize the biological processes
#          underlying the primary ASD-associated module.
#
# Input:  data/processed/ME11_all_genes_annotated.csv
# Output: results/enrichment/ME11_GO_BP.csv
#         results/enrichment/ME11_GO_MF.csv
#         results/enrichment/ME11_KEGG.csv
#         figures/ME11_GO_dotplot.png
#         figures/ME11_KEGG_dotplot.png
#
# Author: [your name]
# Date: 2026-07-17
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# Create output directory
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 1: LOAD ME11 GENES
# ============================================================

cat("Loading ME11 gene list...\n")

ME11_annot <- read.csv("data/processed/ME11_all_genes_annotated.csv")

# Keep characterized genes only
ME11_annot <- ME11_annot[!grepl("^LOC", ME11_annot$Symbol) &
                           ME11_annot$Symbol != "" &
                           !is.na(ME11_annot$Entrez_Gene_ID), ]

# Use Entrez IDs for enrichment (most reliable identifier)
ME11_entrez <- as.character(ME11_annot$Entrez_Gene_ID)
ME11_entrez <- unique(ME11_entrez[!is.na(ME11_entrez) & ME11_entrez != ""])
cat("ME11 genes with Entrez IDs:", length(ME11_entrez), "\n")

# ============================================================
# SECTION 2: GO ENRICHMENT
# ============================================================

cat("\nRunning GO enrichment analysis...\n")

# Biological Process
go_bp <- enrichGO(
  gene          = ME11_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

cat("GO Biological Process terms (FDR<0.05):", nrow(go_bp@result[go_bp@result$p.adjust < 0.05,]), "\n")

# Molecular Function
go_mf <- enrichGO(
  gene          = ME11_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

cat("GO Molecular Function terms (FDR<0.05):", nrow(go_mf@result[go_mf@result$p.adjust < 0.05,]), "\n")

# Save GO results
write.csv(go_bp@result, "results/enrichment/ME11_GO_BP.csv", row.names = FALSE)
write.csv(go_mf@result, "results/enrichment/ME11_GO_MF.csv", row.names = FALSE)
cat("GO results saved.\n")

# Print top terms
cat("\nTop 15 GO Biological Process terms:\n")
print(head(go_bp@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15))

# ============================================================
# SECTION 3: KEGG ENRICHMENT
# ============================================================

cat("\nRunning KEGG pathway enrichment...\n")

kegg <- enrichKEGG(
  gene          = ME11_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20
)

cat("KEGG pathways (FDR<0.05):", nrow(kegg@result[kegg@result$p.adjust < 0.05,]), "\n")

# Convert Entrez IDs to symbols in KEGG results
kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(kegg_readable@result, "results/enrichment/ME11_KEGG.csv", row.names = FALSE)

cat("\nTop KEGG pathways:\n")
print(head(kegg_readable@result[, c("Description", "GeneRatio", "p.adjust", "geneID")], 15))

# ============================================================
# SECTION 4: VISUALIZATIONS
# ============================================================

cat("\nGenerating enrichment plots...\n")

# GO BP dotplot
if(nrow(go_bp@result[go_bp@result$p.adjust < 0.05,]) > 0) {
  p_go <- dotplot(go_bp, showCategory = 20) +
    labs(title = "GO Biological Process — ME11 Module") +
    theme(axis.text.y = element_text(size = 8))
  ggsave("figures/ME11_GO_dotplot.png", p_go, width = 9, height = 8)
  cat("GO dotplot saved.\n")
} else {
  cat("No significant GO BP terms — skipping dotplot.\n")
}

# KEGG dotplot
if(nrow(kegg@result[kegg@result$p.adjust < 0.05,]) > 0) {
  p_kegg <- dotplot(kegg_readable, showCategory = 20) +
    labs(title = "KEGG Pathways — ME11 Module") +
    theme(axis.text.y = element_text(size = 8))
  ggsave("figures/ME11_KEGG_dotplot.png", p_kegg, width = 9, height = 8)
  cat("KEGG dotplot saved.\n")
} else {
  cat("No significant KEGG pathways — skipping dotplot.\n")
}

# GO BP barplot (useful even with few terms)
if(nrow(go_bp) > 0) {
  p_bar <- barplot(go_bp, showCategory = 20) +
    labs(title = "GO Biological Process — ME11 Module") +
    theme(axis.text.y = element_text(size = 8))
  ggsave("figures/ME11_GO_barplot.png", p_bar, width = 9, height = 8)
  cat("GO barplot saved.\n")
}

# ============================================================
# SECTION 5: SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("Enrichment analysis complete.\n")
cat("GO BP terms (FDR<0.05):", nrow(go_bp@result[go_bp@result$p.adjust < 0.05,]), "\n")
cat("GO MF terms (FDR<0.05):", nrow(go_mf@result[go_mf@result$p.adjust < 0.05,]), "\n")
cat("KEGG pathways (FDR<0.05):", nrow(kegg@result[kegg@result$p.adjust < 0.05,]), "\n")
cat("Results saved to results/enrichment/\n")
cat("============================================================\n")

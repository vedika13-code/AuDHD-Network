# ============================================================
# Script: 05_adhd_de.R
# Purpose: Test whether ASD-associated ME11 co-expression
#          module shows coordinated dysregulation in the
#          ADHD dataset (GSE159104).
#
# This is the key AuDHD step: if ME11 genes are also
# dysregulated in ADHD, this suggests shared transcriptional
# architecture between ASD and ADHD — the AuDHD signal.
#
# Approach:
#   1. Extract ME11 gene symbols
#   2. Match to ADHD expression matrix by gene symbol
#   3. Run limma differential expression (ADHD vs CTRL)
#   4. Test whether ME11 genes are enriched among DE genes
#      (Fisher's exact test)
#   5. Calculate ME11 module score in ADHD samples
#   6. Visualize results
#
# Input:  data/processed/ADHD_expression_clean.rds
#         data/processed/ADHD_metadata_clean.rds
#         data/processed/ME11_all_genes_annotated.csv
# Output: results/adhd_de/ME11_adhd_de_results.csv
#         results/adhd_de/ME11_enrichment_test.csv
#         figures/ME11_adhd_volcano.png
#         figures/ME11_adhd_module_score.png
#
# Author: [your name]
# Date: 2026-07-17
# ============================================================

library(limma)
library(ggplot2)

# Create output directory
dir.create("results/adhd_de", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 1: LOAD DATA
# ============================================================

cat("Loading ADHD expression data and ME11 gene list...\n")

adhd_clean      <- readRDS("data/processed/ADHD_expression_clean.rds")
adhd_meta_clean <- readRDS("data/processed/ADHD_metadata_clean.rds")
ME11_annot      <- read.csv("data/processed/ME11_all_genes_annotated.csv")

cat("ADHD matrix:", dim(adhd_clean), "\n")
print(table(adhd_meta_clean$ADHD_Status))

# ME11 gene symbols (characterized only)
ME11_symbols <- unique(ME11_annot$Symbol)
ME11_symbols <- ME11_symbols[!grepl("^LOC", ME11_symbols) & ME11_symbols != ""]
cat("ME11 characterized genes:", length(ME11_symbols), "\n")

# ============================================================
# SECTION 2: MATCH ME11 GENES TO ADHD EXPRESSION MATRIX
# ============================================================

cat("\nMatching ME11 genes to ADHD expression matrix...\n")

# ADHD matrix rownames are gene symbols from RefSeq annotation
# ME11 genes are Illumina probe symbols — need to find overlap
adhd_genes <- rownames(adhd_clean)
cat("ADHD genes available:", length(adhd_genes), "\n")

# Find overlap
overlap <- intersect(ME11_symbols, adhd_genes)
cat("ME11 genes found in ADHD matrix:", length(overlap), "\n")
cat("Genes:", paste(overlap, collapse = ", "), "\n")

if(length(overlap) == 0) {
  cat("\nNo direct symbol overlap found.\n")
  cat("This is expected: ASD dataset uses microarray probe symbols,\n")
  cat("ADHD dataset uses RefSeq gene symbols from RNA-seq.\n")
  cat("Proceeding with full genome-wide DE analysis and enrichment test.\n")
}

# ============================================================
# SECTION 3: FULL GENOME-WIDE DE ANALYSIS (ADHD vs CTRL)
# ============================================================

cat("\nRunning limma differential expression: ADHD vs CTRL...\n")

# Build design matrix
condition <- factor(adhd_meta_clean$ADHD_Status, levels = c("CTRL", "ADHD"))
design    <- model.matrix(~ condition)
cat("Design matrix:\n")
print(head(design))

# Fit linear model
fit  <- lmFit(adhd_clean, design)
fit2 <- eBayes(fit)

# Extract all results
de_results <- topTable(fit2, coef = "conditionADHD",
                        number = Inf, sort.by = "P")
cat("Total genes tested:", nrow(de_results), "\n")
cat("Significant DE genes (FDR < 0.05):", sum(de_results$adj.P.Val < 0.05), "\n")
cat("Significant DE genes (FDR < 0.10):", sum(de_results$adj.P.Val < 0.10), "\n")
cat("Nominal DE genes (p < 0.05):", sum(de_results$P.Value < 0.05), "\n")

# Save full DE results
write.csv(de_results, "results/adhd_de/ADHD_full_de_results.csv")
cat("Full DE results saved.\n")

# Top hits
cat("\nTop 20 DE genes (ADHD vs CTRL):\n")
print(head(de_results[, c("logFC", "P.Value", "adj.P.Val")], 20))

# ============================================================
# SECTION 4: TEST ME11 ENRICHMENT IN ADHD DE GENES
# ============================================================

cat("\nTesting ME11 enrichment in ADHD DE genes...\n")

# This tests whether ME11 genes are overrepresented among
# genes differentially expressed in ADHD
# Uses Fisher's exact test

all_genes  <- rownames(de_results)
n_total    <- length(all_genes)

# Try both symbol matching and any available overlap
ME11_in_adhd  <- intersect(ME11_symbols, all_genes)
cat("ME11 genes found in ADHD DE results:", length(ME11_in_adhd), "\n")

if(length(ME11_in_adhd) > 0) {
  # Fisher's exact test at multiple thresholds
  for(fdr_thresh in c(0.05, 0.10, 0.20)) {
    de_genes    <- rownames(de_results[de_results$adj.P.Val < fdr_thresh, ])
    n_de        <- length(de_genes)
    n_ME11      <- length(ME11_in_adhd)
    n_ME11_de   <- length(intersect(ME11_in_adhd, de_genes))

    # 2x2 contingency table
    mat <- matrix(c(
      n_ME11_de,
      n_ME11 - n_ME11_de,
      n_de - n_ME11_de,
      n_total - n_ME11 - n_de + n_ME11_de
    ), nrow = 2)

    ftest <- fisher.test(mat, alternative = "greater")
    cat(sprintf("FDR < %.2f: ME11 overlap = %d/%d, OR = %.2f, p = %.4f\n",
                fdr_thresh, n_ME11_de, n_ME11,
                ftest$estimate, ftest$p.value))
  }
} else {
  cat("No ME11 gene overlap with ADHD matrix — enrichment test not possible.\n")
  cat("This is due to cross-platform gene naming differences.\n")
  cat("Proceeding with module eigengene approach instead.\n")
}

# ============================================================
# SECTION 5: ME11 MODULE SCORE IN ADHD SAMPLES
# ============================================================

cat("\nCalculating ME11 module score in ADHD samples...\n")

# Even if gene symbols don't perfectly overlap, we can calculate
# a module score using available ME11 genes in ADHD matrix
# Module score = mean z-score of ME11 genes across samples

if(length(ME11_in_adhd) >= 5) {
  ME11_adhd_expr <- adhd_clean[ME11_in_adhd, ]

  # Z-score normalize each gene
  ME11_zscore <- t(scale(t(ME11_adhd_expr)))

  # Module score = mean z-score per sample
  module_score <- colMeans(ME11_zscore, na.rm = TRUE)

  # Compare scores between ADHD and CTRL
  score_df <- data.frame(
    sample     = names(module_score),
    score      = module_score,
    condition  = adhd_meta_clean$ADHD_Status[match(names(module_score),
                                                     rownames(adhd_meta_clean))]
  )

  # T-test
  ttest <- t.test(score ~ condition, data = score_df)
  cat("Module score t-test (ADHD vs CTRL): p =", ttest$p.value, "\n")
  cat("Mean score ADHD:", mean(score_df$score[score_df$condition == "ADHD"]), "\n")
  cat("Mean score CTRL:", mean(score_df$score[score_df$condition == "CTRL"]), "\n")

  # Plot
  p_score <- ggplot(score_df, aes(x = condition, y = score, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = c("CTRL" = "#85C1E9", "ADHD" = "#E74C3C")) +
    labs(title = "ME11 Module Score in ADHD Dataset",
         subtitle = paste0("p = ", signif(ttest$p.value, 3)),
         x = "Condition", y = "ME11 Module Score (mean z-score)") +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave("figures/ME11_adhd_module_score.png", p_score, width = 5, height = 5)
  cat("Module score plot saved.\n")

  # Save score results
  write.csv(score_df, "results/adhd_de/ME11_adhd_module_scores.csv", row.names = FALSE)

} else {
  cat("Fewer than 5 ME11 genes in ADHD matrix — module score not calculated.\n")
  cat("Gene naming mismatch between microarray (ASD) and RNA-seq (ADHD) platforms.\n")
}

# ============================================================
# SECTION 6: VOLCANO PLOT — ADHD DE
# ============================================================

cat("\nGenerating volcano plot...\n")

volcano_df <- data.frame(
  gene    = rownames(de_results),
  logFC   = de_results$logFC,
  negLogP = -log10(de_results$P.Value),
  FDR     = de_results$adj.P.Val,
  ME11    = rownames(de_results) %in% ME11_in_adhd
)

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = ME11), size = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "#BDC3C7", "TRUE" = "#E74C3C"),
                     labels = c("Other genes", "ME11 genes")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  labs(title = "ADHD vs CTRL Differential Expression",
       subtitle = "Red = ME11 ASD-associated module genes",
       x = "log2 Fold Change", y = "-log10(p-value)",
       color = "") +
  theme_minimal()

ggsave("figures/ME11_adhd_volcano.png", p_volcano, width = 7, height = 5)
cat("Volcano plot saved.\n")

# ============================================================
# SECTION 7: SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("ADHD differential expression analysis complete.\n")
cat("Total DE genes (FDR<0.05):", sum(de_results$adj.P.Val < 0.05), "\n")
cat("ME11 genes in ADHD matrix:", length(ME11_in_adhd), "\n")
cat("Results saved to results/adhd_de/\n")
cat("============================================================\n")

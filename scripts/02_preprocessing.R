# ============================================================
# Script: 02_preprocessing.R
# Purpose: Load, clean, normalize, and prepare expression
#          matrices for ADHD (GSE159104) and ASD
#          (GSE111175 + GSE42133 merged) datasets for WGCNA
#
# Key steps:
#   1. Load expression data from GEO
#   2. Subset ASD to ASD + TD only
#   3. Merge ASD datasets + ComBat batch correction
#   4. Filter low-expression genes
#   5. Remove outlier samples
#   6. Rank-based inverse normal transformation (both datasets)
#   7. Save clean matrices
#
# Author: [your name]
# Date: 2026-07-12
# ============================================================

library(GEOquery)
library(limma)
library(sva)
library(WGCNA)
library(ggplot2)

# Allow multi-threading for WGCNA
allowWGCNAThreads()

# ============================================================
# SECTION 1: LOAD EXPRESSION DATA
# ============================================================

cat("Loading expression data from GEO...\n")

# --- ADHD: GSE159104 ---
# Note: This dataset stores counts in supplementary files, not
# the series matrix. The series matrix contains metadata only
# (data_row_count = 0). We load the series matrix for metadata
# and will need to download count data separately.

gse_adhd <- getGEO("GSE159104", GSEMatrix = TRUE)
adhd_meta <- pData(gse_adhd[[1]])

# Download supplementary count file from GEO
# This is the DGE counts file deposited by the authors
adhd_supp <- getGEOSuppFiles("GSE159104", baseDir = "data/raw/")
cat("ADHD supplementary files downloaded:\n")
print(rownames(adhd_supp))

# Load the count matrix (tab-separated)
adhd_count_file <- rownames(adhd_supp)[grepl(".txt|.tsv|.csv", rownames(adhd_supp))]

if(length(adhd_count_file) == 0) {
  # Try gz file
  adhd_count_file <- rownames(adhd_supp)[grepl(".gz", rownames(adhd_supp))]
}

cat("Loading ADHD count file:", adhd_count_file, "\n")
adhd_counts <- read.table(gzfile(adhd_count_file),
                          header = TRUE,
                          row.names = 1,
                          sep = "\t",
                          check.names = FALSE)

cat("ADHD expression matrix dimensions:", dim(adhd_counts), "\n")

# --- ASD: GSE111175 ---
gse_asd1 <- getGEO("GSE111175", GSEMatrix = TRUE)
asd1_expr <- exprs(gse_asd1[[1]])   # expression matrix
asd1_meta <- pData(gse_asd1[[1]])   # metadata

cat("GSE111175 expression matrix dimensions:", dim(asd1_expr), "\n")

# --- ASD: GSE42133 ---
gse_asd2 <- getGEO("GSE42133", GSEMatrix = TRUE)
asd2_expr <- exprs(gse_asd2[[1]])
asd2_meta <- pData(gse_asd2[[1]])

cat("GSE42133 expression matrix dimensions:", dim(asd2_expr), "\n")

# ============================================================
# SECTION 2: SUBSET ASD DATASETS TO ASD + TD ONLY
# ============================================================

cat("\nSubsetting ASD datasets...\n")

# GSE111175: keep ASD and TD only
asd1_keep <- asd1_meta$`diagnosis:ch1` %in% c("ASD", "TD")
asd1_expr  <- asd1_expr[, asd1_keep]
asd1_meta  <- asd1_meta[asd1_keep, ]
cat("GSE111175 after subsetting:", ncol(asd1_expr), "samples\n")
print(table(asd1_meta$`diagnosis:ch1`))

# GSE42133: keep ASD and Control only
asd2_keep <- asd2_meta$`dx (diagnosis):ch1` %in% c("ASD", "Control")
asd2_expr  <- asd2_expr[, asd2_keep]
asd2_meta  <- asd2_meta[asd2_keep, ]
cat("GSE42133 after subsetting:", ncol(asd2_expr), "samples\n")
print(table(asd2_meta$`dx (diagnosis):ch1`))

# Harmonize diagnosis labels in GSE42133 to match GSE111175
# "Control" -> "TD" for consistency
asd2_meta$diagnosis <- ifelse(asd2_meta$`dx (diagnosis):ch1` == "ASD", "ASD", "TD")
asd1_meta$diagnosis <- asd1_meta$`diagnosis:ch1`

# ============================================================
# SECTION 3: MERGE ASD DATASETS
# ============================================================

cat("\nMerging ASD datasets...\n")

# Find common probes between both datasets
common_probes <- intersect(rownames(asd1_expr), rownames(asd2_expr))
cat("Common probes:", length(common_probes), "\n")

# Subset both to common probes
asd1_expr <- asd1_expr[common_probes, ]
asd2_expr <- asd2_expr[common_probes, ]

# Merge expression matrices
asd_merged_expr <- cbind(asd1_expr, asd2_expr)

# Merge metadata
asd_merged_meta <- data.frame(
  sample_id  = c(rownames(asd1_meta), rownames(asd2_meta)),
  diagnosis  = c(asd1_meta$diagnosis, asd2_meta$diagnosis),
  gender     = c(asd1_meta$`gender:ch1`, asd2_meta$`gender:ch1`),
  batch      = c(rep("GSE111175", ncol(asd1_expr)),
                 rep("GSE42133",  ncol(asd2_expr))),
  stringsAsFactors = FALSE
)
rownames(asd_merged_meta) <- asd_merged_meta$sample_id

cat("Merged ASD matrix dimensions:", dim(asd_merged_expr), "\n")
print(table(asd_merged_meta$diagnosis))
print(table(asd_merged_meta$batch))

# ============================================================
# SECTION 4: BATCH CORRECTION (ASD only)
# ============================================================

cat("\nApplying ComBat batch correction to merged ASD data...\n")

# PCA before batch correction
pca_before <- prcomp(t(asd_merged_expr), scale. = TRUE)
pca_before_df <- data.frame(
  PC1   = pca_before$x[, 1],
  PC2   = pca_before$x[, 2],
  batch = asd_merged_meta$batch
)

p_before <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2) +
  labs(title = "PCA Before ComBat — ASD Merged") +
  theme_minimal()

ggsave("figures/pca_asd_before_combat.png", p_before, width = 6, height = 5)

# Apply ComBat
batch_vector <- asd_merged_meta$batch
mod <- model.matrix(~ diagnosis, data = asd_merged_meta)

asd_corrected_expr <- ComBat(
  dat    = asd_merged_expr,
  batch  = batch_vector,
  mod    = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

# PCA after batch correction
pca_after <- prcomp(t(asd_corrected_expr), scale. = TRUE)
pca_after_df <- data.frame(
  PC1   = pca_after$x[, 1],
  PC2   = pca_after$x[, 2],
  batch = asd_merged_meta$batch
)

p_after <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2) +
  labs(title = "PCA After ComBat — ASD Merged") +
  theme_minimal()

ggsave("figures/pca_asd_after_combat.png", p_after, width = 6, height = 5)
cat("PCA plots saved to figures/\n")

# ============================================================
# SECTION 5: PROCESS ADHD EXPRESSION DATA
# ============================================================

cat("\nProcessing ADHD expression data...\n")

# Align ADHD metadata with count matrix columns
# The count matrix columns should match GSM sample IDs
# Check alignment
cat("ADHD count matrix columns (first 5):", colnames(adhd_counts)[1:5], "\n")
cat("ADHD metadata rownames (first 5):", rownames(adhd_meta)[1:5], "\n")

# Keep only ADHD and CTRL samples
adhd_keep <- adhd_meta$`condition:ch1` %in% c("ADHD", "CTRL")
adhd_meta  <- adhd_meta[adhd_keep, ]

# Subset count matrix to matching samples
shared_adhd_samples <- intersect(colnames(adhd_counts), rownames(adhd_meta))
adhd_counts <- adhd_counts[, shared_adhd_samples]
adhd_meta   <- adhd_meta[shared_adhd_samples, ]

cat("ADHD samples after alignment:", ncol(adhd_counts), "\n")
print(table(adhd_meta$`condition:ch1`))

# ============================================================
# SECTION 6: FILTER LOW-EXPRESSION GENES
# ============================================================

cat("\nFiltering low-expression genes...\n")

# --- ASD (microarray): remove probes with low variance ---
# For microarray, filter by variance across samples
# Keep top 75% most variable probes
asd_vars   <- apply(asd_corrected_expr, 1, var)
asd_cutoff <- quantile(asd_vars, 0.25)
asd_filtered_expr <- asd_corrected_expr[asd_vars > asd_cutoff, ]
cat("ASD probes after variance filter:", nrow(asd_filtered_expr), "\n")

# --- ADHD (RNA-seq): filter by expression level ---
# Keep genes expressed in at least 50% of samples with RPKM > 0.5
adhd_counts_num <- apply(adhd_counts, 2, as.numeric)
rownames(adhd_counts_num) <- rownames(adhd_counts)

expressed <- rowSums(adhd_counts_num > 0.5) >= (0.5 * ncol(adhd_counts_num))
adhd_filtered_expr <- adhd_counts_num[expressed, ]
cat("ADHD genes after expression filter:", nrow(adhd_filtered_expr), "\n")

# ============================================================
# SECTION 7: RANK-BASED INVERSE NORMAL TRANSFORMATION
# ============================================================

cat("\nApplying rank-based inverse normal transformation...\n")

# This makes both datasets unit-normal and removes
# platform-specific scaling — critical for cross-platform
# WGCNA module preservation comparison

rankNormTransform <- function(expr_matrix) {
  # Apply rank-based inverse normal transformation per gene
  t(apply(expr_matrix, 1, function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  }))
}

asd_rint  <- rankNormTransform(asd_filtered_expr)
adhd_rint <- rankNormTransform(adhd_filtered_expr)

cat("ASD matrix after RINT:", dim(asd_rint), "\n")
cat("ADHD matrix after RINT:", dim(adhd_rint), "\n")

# ============================================================
# SECTION 8: OUTLIER SAMPLE REMOVAL
# ============================================================

cat("\nRemoving outlier samples...\n")

removeOutliers <- function(expr_matrix, label = "dataset", z_threshold = 2.5) {
  # Hierarchical clustering on sample-sample distances
  sample_tree <- hclust(dist(t(expr_matrix)), method = "average")

  # Calculate standardized connectivity
  connectivity <- softConnectivity(t(expr_matrix), power = 6)
  z_connectivity <- (connectivity - mean(connectivity)) / sd(connectivity)

  # Flag outliers
  outliers <- names(z_connectivity[z_connectivity < -z_threshold])
  cat(label, "- Outlier samples removed:", length(outliers), "\n")
  if(length(outliers) > 0) cat("  Samples:", outliers, "\n")

  # Save clustering plot
  png(paste0("figures/sample_clustering_", label, ".png"),
      width = 1200, height = 600)
  plot(sample_tree, main = paste("Sample Clustering —", label),
       sub = "", xlab = "")
  abline(h = quantile(sample_tree$height, 0.95), col = "red", lty = 2)
  dev.off()

  # Return cleaned matrix
  keep <- setdiff(colnames(expr_matrix), outliers)
  return(expr_matrix[, keep])
}

asd_clean  <- removeOutliers(asd_rint,  label = "ASD",  z_threshold = 2.5)
adhd_clean <- removeOutliers(adhd_rint, label = "ADHD", z_threshold = 2.5)

# Update metadata to match cleaned matrices
asd_meta_clean  <- asd_merged_meta[colnames(asd_clean), ]
adhd_meta_clean <- adhd_meta[colnames(adhd_clean), ]

cat("\nFinal sample counts:\n")
cat("ASD:", ncol(asd_clean), "samples\n")
print(table(asd_meta_clean$diagnosis))
cat("ADHD:", ncol(adhd_clean), "samples\n")
print(table(adhd_meta_clean$`condition:ch1`))

# ============================================================
# SECTION 9: SAVE CLEANED MATRICES
# ============================================================

cat("\nSaving cleaned expression matrices...\n")

# Expression matrices (genes/probes x samples)
write.csv(asd_clean,
          "data/processed/ASD_expression_clean.csv")
write.csv(adhd_clean,
          "data/processed/ADHD_expression_clean.csv")

# Metadata
write.csv(asd_meta_clean,
          "data/processed/ASD_metadata_clean.csv")
write.csv(adhd_meta_clean,
          "data/processed/ADHD_metadata_clean.csv")

# Also save as R objects for faster loading in subsequent scripts
saveRDS(asd_clean,       "data/processed/ASD_expression_clean.rds")
saveRDS(adhd_clean,      "data/processed/ADHD_expression_clean.rds")
saveRDS(asd_meta_clean,  "data/processed/ASD_metadata_clean.rds")
saveRDS(adhd_meta_clean, "data/processed/ADHD_metadata_clean.rds")

cat("\n============================================================\n")
cat("Preprocessing complete.\n")
cat("ASD clean matrix:", nrow(asd_clean), "probes x", ncol(asd_clean), "samples\n")
cat("ADHD clean matrix:", nrow(adhd_clean), "genes x", ncol(adhd_clean), "samples\n")
cat("Output files saved to data/processed/\n")
cat("Figures saved to figures/\n")
cat("============================================================\n")






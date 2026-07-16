# ============================================================
# Script: 02_preprocessing.R
# Purpose: Load, clean, and prepare expression matrices for
#          ASD (GSE111175 + GSE42133 merged) and ADHD
#          (GSE159104) datasets for downstream WGCNA and
#          differential expression analysis
#
# Key decisions documented here:
#   - RINT NOT applied to ASD: destroys microarray correlation
#     structure needed for WGCNA (post-RINT r ~ 0, unusable)
#   - ComBat batch correction applied to ASD merged cohort
#   - ADHD replicates averaged; protein-coding genes only
#   - Outlier removal via sample-sample Pearson correlation
#
# Author: [your name]
# Date: 2026-07-16
# ============================================================

library(GEOquery)
library(sva)
library(WGCNA)
library(ggplot2)

allowWGCNAThreads()

# ============================================================
# SECTION 1: ASD DATASETS (GSE111175 + GSE42133)
# ============================================================

cat("=== Loading ASD datasets ===\n")

gse_asd1 <- getGEO("GSE111175", GSEMatrix = TRUE)
asd1_expr <- exprs(gse_asd1[[1]])
asd1_meta <- pData(gse_asd1[[1]])
cat("GSE111175 loaded:", dim(asd1_expr), "\n")

gse_asd2 <- getGEO("GSE42133", GSEMatrix = TRUE)
asd2_expr <- exprs(gse_asd2[[1]])
asd2_meta <- pData(gse_asd2[[1]])
cat("GSE42133 loaded:", dim(asd2_expr), "\n")

# --- Subset to ASD + TD only ---
asd1_keep <- asd1_meta$`diagnosis:ch1` %in% c("ASD", "TD")
asd1_expr  <- asd1_expr[, asd1_keep]
asd1_meta  <- asd1_meta[asd1_keep, ]
cat("GSE111175 after subset:", dim(asd1_expr), "\n")
print(table(asd1_meta$`diagnosis:ch1`))

asd2_keep <- asd2_meta$`dx (diagnosis):ch1` %in% c("ASD", "Control")
asd2_expr  <- asd2_expr[, asd2_keep]
asd2_meta  <- asd2_meta[asd2_keep, ]
cat("GSE42133 after subset:", dim(asd2_expr), "\n")
print(table(asd2_meta$`dx (diagnosis):ch1`))

# --- Harmonize diagnosis labels ---
asd1_meta$diagnosis <- asd1_meta$`diagnosis:ch1`
asd2_meta$diagnosis <- ifelse(asd2_meta$`dx (diagnosis):ch1` == "ASD", "ASD", "TD")

# --- Merge on common probes ---
common_probes  <- intersect(rownames(asd1_expr), rownames(asd2_expr))
cat("Common probes:", length(common_probes), "\n")

asd1_expr      <- asd1_expr[common_probes, ]
asd2_expr      <- asd2_expr[common_probes, ]
asd_merged     <- cbind(asd1_expr, asd2_expr)

asd_merged_meta <- data.frame(
  sample_id = c(rownames(asd1_meta), rownames(asd2_meta)),
  diagnosis = c(asd1_meta$diagnosis, asd2_meta$diagnosis),
  gender    = c(asd1_meta$`gender:ch1`, asd2_meta$`gender:ch1`),
  batch     = c(rep("GSE111175", ncol(asd1_expr)),
                rep("GSE42133",  ncol(asd2_expr))),
  stringsAsFactors = FALSE
)
rownames(asd_merged_meta) <- asd_merged_meta$sample_id
cat("Merged ASD matrix:", dim(asd_merged), "\n")
print(table(asd_merged_meta$diagnosis))

# --- PCA before ComBat ---
pca_before    <- prcomp(t(asd_merged), scale. = TRUE)
df_before     <- data.frame(PC1 = pca_before$x[,1],
                             PC2 = pca_before$x[,2],
                             batch = asd_merged_meta$batch)
p_before <- ggplot(df_before, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2) +
  labs(title = "PCA Before ComBat — ASD") +
  theme_minimal()
ggsave("figures/pca_asd_before_combat.png", p_before, width = 6, height = 5)

# --- ComBat batch correction ---
cat("Applying ComBat...\n")
mod            <- model.matrix(~ diagnosis, data = asd_merged_meta)
asd_corrected  <- ComBat(dat = asd_merged, batch = asd_merged_meta$batch,
                          mod = mod, par.prior = TRUE, prior.plots = FALSE)

# --- PCA after ComBat ---
pca_after  <- prcomp(t(asd_corrected), scale. = TRUE)
df_after   <- data.frame(PC1 = pca_after$x[,1],
                          PC2 = pca_after$x[,2],
                          batch = asd_merged_meta$batch)
p_after <- ggplot(df_after, aes(PC1, PC2, color = batch)) +
  geom_point(size = 2) +
  labs(title = "PCA After ComBat — ASD") +
  theme_minimal()
ggsave("figures/pca_asd_after_combat.png", p_after, width = 6, height = 5)
cat("PCA plots saved.\n")

# --- Filter top 15,000 most variable probes ---
# Note: RINT NOT applied - destroys correlation structure for WGCNA
asd_vars      <- apply(asd_corrected, 1, var)
top15k        <- names(sort(asd_vars, decreasing = TRUE)[1:15000])
asd_top       <- asd_corrected[top15k, ]
cat("ASD after variance filter:", dim(asd_top), "\n")

# --- Outlier removal ---
# Using sample-sample Pearson correlation Z-scores
# Note: WGCNA softConnectivity not used because it runs on
# the larger dimension (genes) when genes >> samples
sample_cor    <- cor(asd_top, method = "pearson")
mean_cor      <- rowMeans(sample_cor)
z_scores      <- (mean_cor - mean(mean_cor)) / sd(mean_cor)
outliers      <- names(z_scores[z_scores < -2.5])
cat("ASD outliers removed:", length(outliers), "\n")
if(length(outliers) > 0) cat("Samples:", paste(outliers, collapse = ", "), "\n")

png("figures/outlier_ASD.png", width = 900, height = 500)
plot(z_scores, main = "Sample Connectivity — ASD",
     ylab = "Z score", xlab = "Sample index", pch = 19,
     col = ifelse(z_scores < -2.5, "red", "black"))
abline(h = -2.5, col = "red", lty = 2)
dev.off()

asd_clean      <- asd_top[, !colnames(asd_top) %in% outliers]
asd_meta_clean <- asd_merged_meta[colnames(asd_clean), ]

cat("Final ASD matrix:", dim(asd_clean), "\n")
print(table(asd_meta_clean$diagnosis))

# ============================================================
# SECTION 2: ADHD DATASET (GSE159104)
# ============================================================

cat("\n=== Loading ADHD dataset ===\n")

gse_adhd  <- getGEO("GSE159104", GSEMatrix = TRUE)
adhd_meta_raw <- pData(gse_adhd[[1]])

# Download supplementary files
adhd_supp <- getGEOSuppFiles("GSE159104", baseDir = "data/raw/")

# Load count matrix (DGE raw counts)
adhd_count_file <- rownames(adhd_supp)[grepl("DGE_raw", rownames(adhd_supp))]
cat("Loading:", adhd_count_file, "\n")
adhd_counts <- read.table(gzfile(adhd_count_file),
                           header = TRUE, row.names = 1,
                           sep = "\t", check.names = FALSE)
cat("Raw count matrix:", dim(adhd_counts), "\n")

# Load sample mapping file (Channel ID to sample ID)
samples_file <- "data/raw//GSE159104/GSE159104_samples.txt.gz"
samples_map  <- read.table(gzfile(samples_file),
                            header = TRUE, sep = "\t", check.names = FALSE)

# Separate annotation columns (1-7) from expression columns (8+)
adhd_annot <- adhd_counts[, 1:7]
adhd_expr  <- adhd_counts[, 8:ncol(adhd_counts)]
adhd_expr  <- apply(adhd_expr, 2, as.numeric)
rownames(adhd_expr) <- rownames(adhd_counts)

# Filter to protein-coding genes (mRNA entries in RefSeq description)
coding         <- grepl("mRNA", adhd_annot$Description)
adhd_expr      <- adhd_expr[coding, ]
adhd_annot     <- adhd_annot[coding, ]
cat("Protein-coding genes:", nrow(adhd_expr), "\n")

# Match channel IDs to expression matrix
matched        <- intersect(colnames(adhd_expr), samples_map$Channel)
adhd_expr      <- adhd_expr[, matched]
samples_map    <- samples_map[match(matched, samples_map$Channel), ]
cat("Matched channels:", length(matched), "\n")

# Average technical replicates per individual
samples_map$individual <- gsub("_rep[0-9]+$", "", samples_map$`Sample Name`)
individuals    <- unique(samples_map$individual)

adhd_expr_avg  <- sapply(individuals, function(ind) {
  cols <- samples_map$Channel[samples_map$individual == ind]
  if(length(cols) > 1) rowMeans(adhd_expr[, cols]) else adhd_expr[, cols]
})
cat("After averaging replicates:", dim(adhd_expr_avg), "\n")

# Build clean ADHD metadata
adhd_meta_clean <- samples_map[!duplicated(samples_map$individual), ]
rownames(adhd_meta_clean) <- adhd_meta_clean$individual
adhd_meta_clean <- adhd_meta_clean[individuals, ]
cat("ADHD sample counts:\n")
print(table(adhd_meta_clean$ADHD_Status))

# Filter low-expression genes (RPKM > 0.5 in >= 50% of samples)
expressed      <- rowSums(adhd_expr_avg > 0.5) >= (0.5 * ncol(adhd_expr_avg))
adhd_filtered  <- adhd_expr_avg[expressed, ]
cat("ADHD genes after expression filter:", nrow(adhd_filtered), "\n")

# Outlier removal (handle NAs from RPKM data)
sample_cor_adhd <- cor(adhd_filtered, method = "pearson",
                        use = "pairwise.complete.obs")
mean_cor_adhd   <- rowMeans(sample_cor_adhd, na.rm = TRUE)
z_adhd          <- (mean_cor_adhd - mean(mean_cor_adhd)) / sd(mean_cor_adhd)
outliers_adhd   <- names(z_adhd[z_adhd < -2.5])
cat("ADHD outliers removed:", length(outliers_adhd), "\n")

png("figures/outlier_ADHD.png", width = 900, height = 500)
plot(z_adhd, main = "Sample Connectivity — ADHD",
     ylab = "Z score", xlab = "Sample index", pch = 19,
     col = ifelse(z_adhd < -2.5, "red", "black"))
abline(h = -2.5, col = "red", lty = 2)
dev.off()

adhd_clean <- adhd_filtered[, !colnames(adhd_filtered) %in% outliers_adhd]
adhd_meta_clean <- adhd_meta_clean[colnames(adhd_clean), ]

cat("Final ADHD matrix:", dim(adhd_clean), "\n")
print(table(adhd_meta_clean$ADHD_Status))

# ============================================================
# SECTION 3: SAVE OUTPUTS
# ============================================================

cat("\n=== Saving outputs ===\n")

saveRDS(asd_clean,       "data/processed/ASD_expression_clean.rds")
saveRDS(asd_meta_clean,  "data/processed/ASD_metadata_clean.rds")
saveRDS(adhd_clean,      "data/processed/ADHD_expression_clean.rds")
saveRDS(adhd_meta_clean, "data/processed/ADHD_metadata_clean.rds")

# Save metadata as CSV (small enough for git)
write.csv(asd_meta_clean,  "data/processed/ASD_metadata_clean.csv",  row.names = FALSE)
write.csv(adhd_meta_clean, "data/processed/ADHD_metadata_clean.csv", row.names = FALSE)

cat("============================================================\n")
cat("Preprocessing complete.\n")
cat("ASD:  ", nrow(asd_clean),  "probes x", ncol(asd_clean),  "samples\n")
cat("ADHD: ", nrow(adhd_clean), "genes  x", ncol(adhd_clean), "samples\n")
cat("Outputs saved to data/processed/\n")
cat("============================================================\n")

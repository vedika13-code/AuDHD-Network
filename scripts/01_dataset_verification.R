# ============================================================
# Script: 01_dataset_verification.R
# Purpose: Pull metadata for candidate ASD and ADHD datasets
#          from GEO to verify sample size, tissue, platform
# Author: [your name]
# Date: 2026-07-12
# ============================================================

library(GEOquery)

# --- ADHD Dataset: GSE159104 ---
# Whole blood RNA-seq, 143 ADHD + 169 controls
# Sudre et al. (2020), Mol Psychiatry

gse_adhd <- getGEO("GSE159104", GSEMatrix = TRUE)
adhd_meta <- pData(gse_adhd[[1]])

# Preview key columns
print(adhd_meta[, c("title", "characteristics_ch1", "characteristics_ch1.1")])

# Save metadata to file
write.csv(adhd_meta, "data/processed/GSE159104_metadata.csv", row.names = FALSE)
cat("ADHD metadata saved.\n")

# --- ASD Dataset: GSE111175 ---
# Whole blood RNA-seq, ASD vs neurotypical controls

gse_asd <- getGEO("GSE111175", GSEMatrix = TRUE)
asd_meta <- pData(gse_asd[[1]])

print(asd_meta[, c("title", "characteristics_ch1", "characteristics_ch1.1")])

write.csv(asd_meta, "data/processed/GSE111175_metadata.csv", row.names = FALSE)
cat("ASD metadata saved.\n")
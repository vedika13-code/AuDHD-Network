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


# See all column names
colnames(asd_meta)

# See full metadata for first 5 samples
asd_meta[1:5, ] |> t()

gse_asd2 <- getGEO("GSE42133", GSEMatrix = TRUE)
asd2_meta <- pData(gse_asd2[[1]])
colnames(asd2_meta)
asd2_meta[1:5,] |> t()

# Check what platform GSE159104 actually is
adhd_meta <- pData(gse_adhd[[1]])
adhd_meta[1:5,] |> t()

# Count samples per group in ADHD dataset
# First collapse replicates to get unique individuals
table(adhd_meta$`condition:ch1`)

# Count samples per group in ASD dataset
table(asd_meta$`diagnosis:ch1`)

table(asd2_meta$`dx:ch1`)

colnames(asd2_meta)

table(asd2_meta$`dx (diagnosis):ch1`)

asd_meta$platform_id[1]
asd2_meta$platform_id[1]

# Check age range in both to make sure populations are comparable
summary(as.numeric(asd_meta$`age (months):ch1`))
asd2_meta[1:3,] |> t()

colnames(asd2_meta)
# look for an age column, then:
summary(as.numeric(asd2_meta$`age:ch1`))  # adjust column name if needed

# Check if GSE159104 has an age column
colnames(adhd_meta)

# Find age in GSE42133 - check all characteristics columns
asd2_meta[1:3, grep("ch1", colnames(asd2_meta))]

# Also check GSE159104 - maybe age is buried in characteristics
adhd_meta[1:3, grep("characteristics", colnames(adhd_meta))]

gse_asd3 <- getGEO("GSE212645", GSEMatrix = TRUE)
asd3_meta <- pData(gse_asd3[[1]])
colnames(asd3_meta)
table(asd3_meta$`diagnosis:ch1`)  # adjust column name after seeing colnames

# Check the characteristics columns directly
asd3_meta[1:3, grep("characteristics", colnames(asd3_meta))] |> t()

# Also check age and diagnosis-looking columns
table(asd3_meta$`genotype:ch1`)
table(asd3_meta$`cell type:ch1`)
asd3_meta$`puberty:ch1`[1:5]
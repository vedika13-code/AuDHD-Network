# ============================================================
# Script: 03_wgcna_asd.R
# Purpose: Construct weighted gene co-expression network for
#          ASD dataset and identify modules associated with
#          ASD diagnosis
#
# Key decisions:
#   - Power = 12 (R²=0.88, first stable plateau)
#   - Signed network (preserves direction of correlation)
#   - ComBat-corrected values used directly (no RINT)
#   - RINT was tested and rejected: destroys correlation
#     structure in microarray data (post-RINT r ~ 0)
#   - Module-trait correlation with binary ASD diagnosis
#   - Hub genes: |MM| > 0.8 AND |GS| > 0.2
#
# Input:  data/processed/ASD_expression_clean.rds
#         data/processed/ASD_metadata_clean.rds
# Output: data/processed/ASD_network_final.rds
#         data/processed/ASD_module_eigengenes.rds
#         data/processed/ASD_module_trait_results.rds
#         data/processed/ME11_hub_genes_annotated.csv
#         data/processed/ME11_all_genes_annotated.csv
#         figures/ASD_soft_threshold.png
#         figures/ASD_dendrogram_final.png
#         figures/ASD_module_trait_heatmap.png
#
# Author: [your name]
# Date: 2026-07-16
# ============================================================

library(WGCNA)
library(GEOquery)

allowWGCNAThreads()

# ============================================================
# SECTION 1: LOAD DATA
# ============================================================

cat("Loading ASD data...\n")

asd_clean      <- readRDS("data/processed/ASD_expression_clean.rds")
asd_meta_clean <- readRDS("data/processed/ASD_metadata_clean.rds")

cat("ASD matrix:", dim(asd_clean), "\n")
print(table(asd_meta_clean$diagnosis))

# Transpose: WGCNA expects samples x genes
asd_wgcna <- t(asd_clean)
cat("Transposed for WGCNA:", dim(asd_wgcna), "\n")

# Check for bad genes/samples
gsg <- goodSamplesGenes(asd_wgcna, verbose = 3)
cat("All OK:", gsg$allOK, "\n")
if (!gsg$allOK) {
  asd_wgcna <- asd_wgcna[gsg$goodSamples, gsg$goodGenes]
  cat("After cleaning:", dim(asd_wgcna), "\n")
}

# ============================================================
# SECTION 2: SOFT-THRESHOLDING POWER SELECTION
# ============================================================

cat("\nRunning soft-thresholding power selection...\n")

powers <- c(1:10, seq(12, 20, by = 2))

sft <- pickSoftThreshold(
  asd_wgcna,
  powerVector = powers,
  verbose     = 5,
  networkType = "signed"
)

# Plot
png("figures/ASD_soft_threshold.png", width = 900, height = 500)
par(mfrow = c(1, 2))

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R²)",
     type = "n", main = "Scale Independence — ASD")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.85, col = "blue", lty = 2)

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity — ASD")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers, col = "red")
dev.off()

cat("Soft threshold plot saved.\n")
cat("Automatic power estimate:", sft$powerEstimate, "\n")
cat("Selected power: 12 (R²=0.88, first stable plateau)\n")

# Power selection rationale:
# R² reaches 0.858 at power 10, 0.88 at power 12, plateaus thereafter
# Power 12 selected: R² well above 0.85, mean connectivity (30.9)
# sufficient for meaningful module detection

selected_power <- 12

# ============================================================
# SECTION 3: NETWORK CONSTRUCTION
# ============================================================

cat("\nConstructing ASD co-expression network (power =", selected_power, ")...\n")
cat("This will take 15-20 minutes...\n")

net_asd <- blockwiseModules(
  asd_wgcna,
  power             = selected_power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = TRUE,
  saveTOMFileBase   = "data/processed/ASD_TOM",
  verbose           = 3
)

cat("\nNetwork construction complete.\n")
cat("Modules detected:", length(table(net_asd$colors)), "\n")
print(table(net_asd$colors))

# Save network
saveRDS(net_asd, "data/processed/ASD_network_final.rds")
cat("Network saved.\n")

# ============================================================
# SECTION 4: VISUALIZE DENDROGRAM
# ============================================================

png("figures/ASD_dendrogram_final.png", width = 1200, height = 600)
plotDendroAndColors(
  net_asd$dendrograms[[1]],
  net_asd$colors[net_asd$blockGenes[[1]]],
  "Module colors",
  dendroLabels  = FALSE,
  hang          = 0.03,
  addGuide      = TRUE,
  guideHang     = 0.05,
  main          = "ASD Gene Dendrogram and Module Colors"
)
dev.off()
cat("Dendrogram saved.\n")

# ============================================================
# SECTION 5: MODULE-TRAIT ASSOCIATION
# ============================================================

cat("\nCalculating module-trait associations...\n")

# Calculate module eigengenes
MEs <- moduleEigengenes(asd_wgcna, net_asd$colors)$eigengenes
MEs <- orderMEs(MEs)

# Binary trait: ASD=1, TD=0
diagnosis_bin <- as.numeric(asd_meta_clean$diagnosis == "ASD")

# Correlate eigengenes with trait
ME_cor  <- cor(MEs, diagnosis_bin, use = "pairwise.complete.obs")
ME_pval <- corPvalueStudent(ME_cor, nrow(asd_wgcna))

# Results table
ME_results <- data.frame(
  Module      = rownames(ME_cor),
  Correlation = ME_cor[,1],
  Pvalue      = ME_pval[,1],
  FDR         = p.adjust(ME_pval[,1], method = "BH")
)
ME_results <- ME_results[order(ME_results$Pvalue), ]

cat("Top modules associated with ASD diagnosis:\n")
print(head(ME_results, 10))

# Save results
write.csv(ME_results, "data/processed/ASD_module_trait.csv", row.names = FALSE)
saveRDS(MEs, "data/processed/ASD_module_eigengenes.rds")
saveRDS(ME_results, "data/processed/ASD_module_trait_results.rds")

# Heatmap
png("figures/ASD_module_trait_heatmap.png", width = 800, height = 1200)
labeledHeatmap(
  Matrix        = ME_cor,
  xLabels       = "ASD Diagnosis",
  yLabels       = rownames(ME_cor),
  ySymbols      = rownames(ME_cor),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = paste0(signif(ME_cor, 2), "\n(p=", signif(ME_pval, 1), ")"),
  setStdMargins = FALSE,
  cex.text      = 0.5,
  zlim          = c(-1, 1),
  main          = "Module-Trait Correlation — ASD"
)
dev.off()
cat("Heatmap saved.\n")

# ============================================================
# SECTION 6: HUB GENE IDENTIFICATION (ME11)
# ============================================================

cat("\nIdentifying hub genes in ME11...\n")

# Module membership and gene significance
gene_MM <- cor(asd_wgcna, MEs, use = "pairwise.complete.obs")
gene_GS <- cor(asd_wgcna, diagnosis_bin, use = "pairwise.complete.obs")

# ME11 genes
ME11_genes <- names(net_asd$colors[net_asd$colors == 11])
cat("Total genes in ME11:", length(ME11_genes), "\n")

MM_ME11 <- gene_MM[ME11_genes, "ME11"]
GS_ME11 <- gene_GS[ME11_genes, 1]

# Hub gene criteria: |MM| > 0.8 AND |GS| > 0.2
hub_genes <- ME11_genes[abs(MM_ME11) > 0.8 & abs(GS_ME11) > 0.2]
cat("Hub genes identified:", length(hub_genes), "\n")

# Annotate using GPL10558
gpl    <- getGEO("GPL10558")
annot  <- Table(gpl)

hub_annot  <- annot[annot$ID %in% hub_genes,
                    c("ID", "Symbol", "ILMN_Gene", "Entrez_Gene_ID")]
ME11_annot <- annot[annot$ID %in% ME11_genes,
                    c("ID", "Symbol", "ILMN_Gene", "Entrez_Gene_ID")]

cat("\nHub genes annotated:\n")
print(hub_annot)

# Save
write.csv(hub_annot,  "data/processed/ME11_hub_genes_annotated.csv",  row.names = FALSE)
write.csv(ME11_annot, "data/processed/ME11_all_genes_annotated.csv",  row.names = FALSE)

cat("\n============================================================\n")
cat("WGCNA complete.\n")
cat("Primary ASD module: ME11 (r=-0.26, FDR=0.0008, n=298 genes)\n")
cat("Hub genes for PPI: TMSB4X, RAP1B, HSPCAL3, C1D, ZCCHC10\n")
cat("============================================================\n")

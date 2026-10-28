# =====================================================================
# 00_make_inputs.R  (PATCH: new file, replaces README step 0)
# Builds the five inputs the revision scripts need from the files the
# existing AuDHD-Network pipeline already saved, so the WGCNA script
# does not have to be edited or re-run.
#
#   inputs/asd_expr_combat.rds        238 x 15000, the matrix WGCNA used
#   inputs/asd_expr_precombat.rds     same samples/probes, before ComBat
#   inputs/asd_expr_precombat_full.rds all common probes, before ComBat (for xCell)
#   inputs/asd_pheno.rds              sample, diagnosis, dataset
#   inputs/module_labels.rds          named numeric vector, 11 = ME11
#   inputs/probe2gene.rds             probe, symbol, entrez (GPL10558)
#
# It stops with an error unless the published numbers are reproduced.
# Run from this folder:  Rscript 00_make_inputs.R
# =====================================================================
project <- path.expand("~/Desktop/AuDHD_Network")   # EDIT if the project folder moves
pp <- function(...) file.path(project, ...)

suppressPackageStartupMessages({ library(GEOquery); library(WGCNA) })
dir.create("inputs", showWarnings = FALSE)

# ---- 1. ComBat matrix, phenotype, module labels (already saved by the pipeline) ----
expr <- readRDS(pp("data/processed/ASD_wgcna_input.rds"))          # samples x probes
meta <- readRDS(pp("data/processed/ASD_metadata_clean.rds"))
net  <- readRDS(pp("data/processed/ASD_network_final.rds"))
meta <- meta[match(rownames(expr), meta$sample_id), ]
labs <- net$colors
if (is.null(names(labs))) names(labs) <- colnames(expr)
labs <- labs[colnames(expr)]
pheno <- data.frame(sample = meta$sample_id, diagnosis = meta$diagnosis, dataset = meta$batch)

# ---- 2. Probe-to-gene table from the cached GPL10558 annotation ----
gpl_file <- pp("data/raw/geo_cache/GPL10558.soft.gz")
gpl <- if (file.exists(gpl_file)) getGEO(filename = gpl_file) else getGEO("GPL10558", destdir = pp("data/raw/geo_cache"))
gt  <- Table(gpl)
p2g <- data.frame(probe = gt$ID, symbol = gt$Symbol, entrez = as.character(gt$Entrez_Gene_ID))

# ---- 3. Pre-ComBat matrix, rebuilt exactly as in 02_preprocessing.R ----
gse <- function(id) {
  f <- pp("data/raw/geo_cache", paste0(id, "_series_matrix.txt.gz"))
  if (file.exists(f)) getGEO(filename = f, getGPL = FALSE)
  else getGEO(id, GSEMatrix = TRUE, getGPL = FALSE, destdir = pp("data/raw/geo_cache"))[[1]]
}
g1 <- gse("GSE111175"); g2 <- gse("GSE42133")
e1 <- exprs(g1); e2 <- exprs(g2)
cp <- intersect(rownames(e1), rownames(e2))
raw <- t(cbind(e1[cp, ], e2[cp, ]))                                   # samples x probes
if (max(raw, na.rm = TRUE) > 100) raw <- log2(pmax(raw, 0) + 1)       # make sure it is log2
raw <- raw[rownames(expr), ]
pre <- raw[, colnames(expr)]

# ---- 4. Checks against the published numbers ----
chk <- function(ok, msg) { if (!isTRUE(ok)) stop("CHECK FAILED: ", msg, call. = FALSE); cat("  ok:", msg, "\n") }
cat("Checks:\n")
chk(identical(dim(expr), c(238L, 15000L)), "ComBat matrix is 238 x 15000")
chk(all(pheno$sample == rownames(expr)), "phenotype rows match the expression samples")
chk(all(table(pheno$diagnosis) == c(115, 123)), "115 ASD and 123 TD")
tb <- table(pheno$dataset, pheno$diagnosis)
chk(tb["GSE111175", "ASD"] == 28 && tb["GSE111175", "TD"] == 69 &&
    tb["GSE42133", "ASD"] == 87 && tb["GSE42133", "TD"] == 54, "28/69 and 87/54 by dataset")
chk(sum(labs == 11) == 298, "ME11 has 298 probes")
chk(identical(dim(pre), dim(expr)) && !anyNA(pre), "pre-ComBat matrix has the same samples and probes, no missing values")
chk(length(cp) == 47323, "47,323 common probes before filtering")
me11 <- moduleEigengenes(expr, labs)$eigengenes[, "ME11"]
r <- cor(me11, as.numeric(pheno$diagnosis == "ASD"))
chk(abs(r - (-0.2603)) < 0.001, sprintf("ME11 vs diagnosis r = %.4f (paper: -0.26)", r))
# ComBat applies one shift/scale per probe per dataset, so within a dataset the
# pre- and post-ComBat values of a probe must be almost perfectly correlated
s <- pheno$dataset == "GSE42133"; pr <- colnames(expr)[1:200]
rr <- sapply(pr, function(p) cor(pre[s, p], expr[s, p]))
chk(median(rr) > 0.99, sprintf("pre-ComBat matches ComBat within dataset (median r = %.4f)", median(rr)))
sym <- p2g$symbol[match(colnames(expr), p2g$probe)]
g11 <- unique(sym[labs == 11 & !is.na(sym) & sym != "" & !grepl("^LOC\\d+", sym)])
chk(length(g11) == 187, sprintf("ME11 maps to %d annotated genes (paper: 187)", length(g11)))

# ---- 5. Save ----
saveRDS(expr,  "inputs/asd_expr_combat.rds")
saveRDS(pre,   "inputs/asd_expr_precombat.rds")
saveRDS(raw,   "inputs/asd_expr_precombat_full.rds")
saveRDS(pheno, "inputs/asd_pheno.rds")
saveRDS(labs,  "inputs/module_labels.rds")
saveRDS(p2g,   "inputs/probe2gene.rds")
cat("All checks passed. Inputs written to inputs/\n")

# =====================================================================
# 00_config.R  -  shared settings and helpers for the revision analyses
# Every other script starts with source("00_config.R").
# Edit the paths in the "EDIT ME" block, then run the scripts in order.
# =====================================================================

# ---------------- EDIT ME ------------------------------------------------
# Objects saved from your existing pipeline (see README, step 0).
f_expr_combat    <- "inputs/asd_expr_combat.rds"     # samples x probes, 238 x 15000, after ComBat (the matrix WGCNA was run on)
f_expr_precombat <- "inputs/asd_expr_precombat.rds"  # same 238 samples and 15000 probes, log2, BEFORE ComBat
f_pheno          <- "inputs/asd_pheno.rds"           # data.frame: sample, diagnosis ("ASD"/"TD"), dataset ("GSE111175"/"GSE42133")
f_labels         <- "inputs/module_labels.rds"       # named vector: names = probe IDs, values = module number (0 = grey, 11 = ME11)
f_probe2gene     <- "inputs/probe2gene.rds"          # data.frame: probe, symbol, entrez (from GPL10558)
f_expr_precombat_full <- "inputs/asd_expr_precombat_full.rds"  # PATCH: all 47,323 probes, before ComBat (used by 03)
# PATCH: all inputs above are written by 00_make_inputs.R (run it first)

target_module    <- 11        # ME11
soft_power       <- 12
out_dir          <- "results"
n_perm_pres      <- 200       # modulePreservation permutations
n_perm_ppi       <- 10000     # PPI null gene sets
# PATCH: 01 now reports two ASD definitions side by side: "broad" (autism + PDD-NOS +
#        Asperger's) and "strict" (autism only, as GSE111175 was handled in the paper).
max_block        <- 8000      # PATCH: genes per block in 02 (20000 = one block, needs ~6 GB RAM)
# --------------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))   # PATCH: so need() can install missing packages
set.seed(2026)

need <- function(pkgs, bioc = FALSE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(miss, update = FALSE, ask = FALSE)
    } else install.packages(miss)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# Collapse a samples x features matrix to one column per gene symbol,
# keeping the feature with the highest mean expression (WGCNA "MaxMean").
collapse_max_mean <- function(X, sym) {
  keep <- !is.na(sym) & sym != "" & !grepl("^LOC\\d+", sym)
  X <- X[, keep, drop = FALSE]; sym <- sym[keep]
  mu <- colMeans(X, na.rm = TRUE)
  o <- order(sym, -mu)
  idx <- o[!duplicated(sym[o])]
  out <- X[, idx, drop = FALSE]
  colnames(out) <- sym[idx]
  attr(out, "source_feature") <- colnames(X)[idx]
  out
}

# Cohen's d (ASD minus TD) with Hedges correction and sampling variance
cohens_d <- function(x, g) {
  a <- x[g == "ASD"]; b <- x[g == "TD"]
  n1 <- sum(!is.na(a)); n2 <- sum(!is.na(b))
  sp <- sqrt(((n1 - 1) * var(a, na.rm = TRUE) + (n2 - 1) * var(b, na.rm = TRUE)) / (n1 + n2 - 2))
  d <- (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g_ <- d * J
  v <- (n1 + n2) / (n1 * n2) + g_^2 / (2 * (n1 + n2))
  data.frame(n_ASD = n1, n_TD = n2, d = g_, var_d = v,
             ci_lo = g_ - 1.96 * sqrt(v), ci_hi = g_ + 1.96 * sqrt(v))
}

# Module score = mean z-score of the module genes present (same as in the paper)
module_score <- function(Xgene, genes) {
  g <- intersect(genes, colnames(Xgene))
  if (length(g) < 10) warning("Only ", length(g), " module genes present")
  rowMeans(scale(Xgene[, g, drop = FALSE]), na.rm = TRUE)
}

# Reference objects used by several scripts
load_reference <- function() {
  expr  <- readRDS(f_expr_combat)
  pheno <- readRDS(f_pheno)
  labs  <- readRDS(f_labels)
  p2g   <- readRDS(f_probe2gene)
  stopifnot(all(rownames(expr) == pheno$sample) || all(rownames(expr) %in% pheno$sample))
  pheno <- pheno[match(rownames(expr), pheno$sample), ]
  labs  <- labs[colnames(expr)]
  sym   <- p2g$symbol[match(colnames(expr), p2g$probe)]
  Xg    <- collapse_max_mean(expr, sym)
  gene_labels <- setNames(labs[attr(Xg, "source_feature")], colnames(Xg))
  # PATCH: module gene sets are defined from ALL probes in each module, as in the paper
  # (ME11 = 187 genes). The one-probe-per-gene collapse above keeps the most highly
  # expressed probe of each gene, and for 87 of the 187 ME11 genes that probe sits in
  # another module, so gene_labels alone would shrink ME11 to 100 genes.
  # gene_labels (one label per gene) is still used where a single label is required
  # (modulePreservation, consensus overlap).
  ok <- !is.na(sym) & sym != "" & !grepl("^LOC\\d+", sym) & labs != 0
  module_sets <- lapply(split(sym[ok], labs[ok]), unique)
  list(expr = expr, pheno = pheno, labels = labs, p2g = p2g,
       Xgene = Xg, gene_labels = gene_labels, module_sets = module_sets,
       module_genes = module_sets[[as.character(target_module)]])
}

# =====================================================================
# 03_cell_composition.R
# Replace marker-average cell scores with xCell deconvolution, run
# separately within each dataset (pre-ComBat), then re-test ME11.
# Also writes CIBERSORTx mixture files if you prefer the LM22 web tool.
# Fixes reviewer point 2 (is ME11 just a T-cell proportion signal?).
# Outputs: results/03_*.csv, results/03_cibersortx_mixture_*.txt
# =====================================================================
source("00_config.R")
need(c("WGCNA", "metafor"))
if (!requireNamespace("xCell", quietly = TRUE)) {
  need("remotes"); remotes::install_github("dviraran/xCell")
}
library(xCell)

ref   <- load_reference()
# PATCH: xCell scores ~10,800 signature genes, so use ALL probes (pre-ComBat),
#        not only the 15,000 most variable ones the network was built on
Xpre  <- readRDS(f_expr_precombat_full)
pheno <- ref$pheno
sym   <- ref$p2g$symbol[match(colnames(Xpre), ref$p2g$probe)]

# ME11 eigengene from the network matrix (as in the paper)
Xm   <- ref$expr[, ref$labels == target_module]
me11 <- moduleEigengenes(Xm, colors = rep(1, ncol(Xm)))$eigengenes[, 1]   # PATCH: colors argument is required
# orient so that it correlates positively with mean module expression
if (cor(me11, rowMeans(scale(ref$expr[, ref$labels == target_module]))) < 0) me11 <- -me11
names(me11) <- rownames(ref$expr)

cells <- c("CD4+ T-cells", "CD8+ T-cells", "NK cells", "B-cells",
           "Monocytes", "Neutrophils", "Eosinophils")
rows <- list(); scores_all <- list()
for (s in unique(pheno$dataset)) {
  idx <- pheno$dataset == s
  Xg  <- collapse_max_mean(Xpre[pheno$sample[idx], , drop = FALSE], sym)
  # CIBERSORTx wants genes x samples on a linear scale
  mix <- data.frame(Gene = colnames(Xg), t(2^Xg), check.names = FALSE)
  write.table(mix, file.path(out_dir, paste0("03_cibersortx_mixture_", s, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  xc  <- t(xCellAnalysis(t(Xg), rnaseq = FALSE, parallel.sz = 2))
  xc  <- xc[, intersect(cells, colnames(xc)), drop = FALSE]
  scores_all[[s]] <- data.frame(sample = rownames(xc), dataset = s, xc, check.names = FALSE)

  df <- data.frame(me = me11[rownames(xc)],
                   dx = factor(pheno$diagnosis[match(rownames(xc), pheno$sample)], c("TD", "ASD")),
                   scale(xc))
  colnames(df) <- make.names(colnames(df))
  f0 <- lm(me ~ dx, df)
  f1 <- lm(as.formula(paste("me ~ dx +", paste(setdiff(colnames(df), c("me", "dx")), collapse = " + "))), df)
  sdy <- sd(df$me)
  for (mod in list(list("unadjusted", f0), list("xCell-adjusted", f1))) {
    co <- summary(mod[[2]])$coefficients["dxASD", ]
    rows[[length(rows) + 1]] <- data.frame(dataset = s, model = mod[[1]],
      beta_SD = co[1] / sdy, se_SD = co[2] / sdy, p = co[4], R2 = summary(mod[[2]])$r.squared)
  }
  # which cell types track ME11?
  cc <- sapply(setdiff(colnames(df), c("me", "dx")), function(v) cor(df$me, df[[v]], method = "spearman"))
  rows[[length(rows) + 1]] <- data.frame(dataset = s, model = paste0("rho_", names(cc)),
                                         beta_SD = cc, se_SD = NA, p = NA, R2 = NA)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "03_ME11_cell_adjusted.csv"), row.names = FALSE)
write.csv(do.call(rbind, scores_all), file.path(out_dir, "03_xcell_scores.csv"), row.names = FALSE)

# meta-analysis of the adjusted effect across the two datasets
adj <- res[res$model == "xCell-adjusted", ]
f <- metafor::rma(yi = adj$beta_SD, sei = adj$se_SD, method = "FE")
write.csv(data.frame(beta_SD = as.numeric(f$b), ci_lo = f$ci.lb, ci_hi = f$ci.ub, p = f$pval, Q_p = f$QEp),
          file.path(out_dir, "03_ME11_cell_adjusted_meta.csv"), row.names = FALSE)
print(res); print(f)

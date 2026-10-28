# =====================================================================
# 01_asd_replication.R
# Does ME11 replicate in independent ASD blood cohorts?
#   (a) ME11 score ASD vs TD in each new cohort (effect size)
#   (b) modulePreservation of every ASD module in each new cohort
#   (c) meta-analysis of ME11 across original + new cohorts
# Fixes reviewer point 2 (single-cohort signal) and part of point 1.
# Outputs: results/01_*.csv
# =====================================================================
source("00_config.R")
need(c("GEOquery", "limma"), bioc = TRUE)
need(c("WGCNA", "metafor"))
allowWGCNAThreads()

options(timeout = 900)                                    # PATCH: allow large GEO files, but give up on a hung download after 15 min
dir.create("geo_cache", showWarnings = FALSE)             # PATCH: cache downloads so a rerun does not start over
ref <- load_reference()
cat("ME11 genes (annotated):", length(ref$module_genes), "\n")

# Candidate replication cohorts. Check each on GEO before relying on it.
cohorts <- c("GSE18123", "GSE25507", "GSE26415")
# PATCH: run only some cohorts with  Rscript 01_asd_replication.R GSE18123
#        results already saved for the other cohorts are kept and the meta-analysis uses all of them
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) cohorts <- args

guess_symbol <- function(fd) {
  for (col in c("Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "gene_symbol")) {
    if (col %in% colnames(fd)) return(trimws(sub(" ///.*", "", fd[[col]])))
  }
  if ("gene_assignment" %in% colnames(fd)) {                 # Affymetrix Gene ST
    return(trimws(vapply(strsplit(fd$gene_assignment, " // "),
                         function(z) if (length(z) >= 2) z[2] else NA_character_, "")))
  }
  stop("No gene-symbol column found; columns are: ", paste(colnames(fd), collapse = ", "))
}

guess_dx <- function(pd) {
  cols <- grep("characteristics|source_name|title|description", colnames(pd), value = TRUE)
  # PATCH: if GEO has a parsed diagnosis/group field, read ONLY that field. In GSE18123 the
  #        word "autism" appears in every sample's description, so every control was called ASD.
  dxcols <- grep("^(diagnosis|dx|group|disease|disease state|sample type|phenotype)[^:]*:ch1$",
                 colnames(pd), ignore.case = TRUE, value = TRUE)
  if (length(dxcols)) cols <- dxcols
  txt <- tolower(apply(pd[, cols, drop = FALSE], 1, paste, collapse = " | "))
  dx <- rep(NA_character_, length(txt))
  dx[grepl("autis|\\basd\\b|pdd|asperger", txt)] <- "ASD"
  # PATCH: \\btypical so that "atypical autism" is not read as a control
  dx[grepl("control|\\btypical|\\btd\\b|healthy|normal", txt) &
     !grepl("autis|\\basd\\b|pdd|asperger", txt)] <- "TD"
  dx[grepl("mother|sibling|parent", txt)] <- NA          # drop relatives
  # PATCH: GSE26415 also contains 42 mothers, labelled "healthy women who had children with ASD (asdMO)"
  #        and "control (ctrlMO)"; without this line they were counted as 21 extra ASD and 21 extra controls
  dx[grepl("women who had child|children with|asdmo|ctrlmo", txt)] <- NA
  # PATCH: strict definition = autism only (PDD-NOS and Asperger's dropped), as in GSE111175
  dx_strict <- dx
  dx_strict[dx == "ASD" & grepl("pdd|asperger", txt)] <- NA
  data.frame(sample = rownames(pd), dx = dx, dx_strict = dx_strict, text = substr(txt, 1, 140))
}

infer_sex <- function(Xg) {
  y <- intersect(c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"), colnames(Xg))
  if (!length(y) || !"XIST" %in% colnames(Xg)) return(rep(NA, nrow(Xg)))
  s <- rowMeans(Xg[, y, drop = FALSE]) - Xg[, "XIST"]
  km <- kmeans(s, 2)
  ifelse(km$cluster == which.max(km$centers), "M", "F")
}

res_effect <- list(); res_pres <- list()

for (gse in cohorts) {
  # PATCH: retry up to 3 times and print the reason when a download fails
  gl <- NULL
  for (attempt in 1:3) {
    gl <- tryCatch(GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE, destdir = "geo_cache"),
                   error = function(e) { message(gse, " attempt ", attempt, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(gl)) break
  }
  if (is.null(gl)) { message("Could not download ", gse); next }
  for (k in seq_along(gl)) {
    es  <- gl[[k]]
    tag <- paste0(gse, "_", Biobase::annotation(es))
    X   <- t(Biobase::exprs(es))                               # samples x probes
    if (max(X, na.rm = TRUE) > 100) X <- log2(pmax(X, 0) + 1)
    sym <- guess_symbol(Biobase::fData(es))
    Xg  <- collapse_max_mean(X, sym)

    dx <- guess_dx(Biobase::pData(es))
    write.csv(dx, file.path(out_dir, paste0("01_check_diagnosis_", tag, ".csv")), row.names = FALSE)
    # >>> Open that CSV and confirm the ASD/TD calls before trusting the results. <<<
    keep <- !is.na(dx$dx)
    Xg <- Xg[keep, , drop = FALSE]; dx <- dx[keep, ]
    if (length(table(dx$dx)) < 2 || min(table(dx$dx)) < 10) { message(tag, ": fewer than 10 per group, skipped"); next }
    sex <- infer_sex(Xg)
    cat(tag, ":", table(dx$dx), " sex:", table(sex), "\n")

    # (a) ME11 score effect, all samples and males only (original cohorts were all male)
    # PATCH: run for both ASD definitions; subset is "all_broad", "males_strict", etc.
    n_genes <- length(intersect(ref$module_genes, colnames(Xg)))
    for (def in c("broad", "strict")) {
      g  <- if (def == "broad") dx$dx else dx$dx_strict
      k  <- !is.na(g)
      if (length(table(g[k])) < 2 || min(table(g[k])) < 10) next
      sc <- module_score(Xg[k, , drop = FALSE], ref$module_genes)   # z-scored within the samples used
      sx <- sex[k]; gk <- g[k]
      df  <- data.frame(sc = sc, dx = factor(gk, levels = c("TD", "ASD")), sex = sx)
      fml <- if (length(unique(na.omit(sx))) > 1) sc ~ dx + sex else sc ~ dx
      p_all <- summary(lm(fml, data = df))$coefficients["dxASD", 4]
      res_effect[[length(res_effect) + 1]] <- cbind(cohort = tag, subset = paste0("all_", def),
                                                    n_ME11_genes = n_genes, cohens_d(sc, gk), p_adj_sex = p_all)
      m <- sx == "M" & !is.na(sx)
      if (sum(m & gk == "ASD") >= 10 && sum(m & gk == "TD") >= 10) {
        res_effect[[length(res_effect) + 1]] <- cbind(cohort = tag, subset = paste0("males_", def),
                                                      n_ME11_genes = n_genes, cohens_d(sc[m], gk[m]),
                                                      p_adj_sex = t.test(sc[m] ~ gk[m])$p.value)
      }
    }

    # (b) preservation of all ASD modules
    common <- intersect(colnames(ref$Xgene), colnames(Xg))
    mp <- modulePreservation(
      multiData  = list(ref  = list(data = ref$Xgene[, common]),
                        test = list(data = Xg[, common])),
      multiColor = list(ref = ref$gene_labels[common]),
      networkType = "signed", referenceNetworks = 1,
      nPermutations = n_perm_pres, randomSeed = 1, quickCor = 0, verbose = 0)
    z  <- mp$preservation$Z$ref.ref$inColumnsAlsoPresentIn.test
    rk <- mp$preservation$observed$ref.ref$inColumnsAlsoPresentIn.test
    res_pres[[tag]] <- data.frame(cohort = tag, module = rownames(z), moduleSize = z$moduleSize,
                                  Zsummary = z$Zsummary.pres, medianRank = rk$medianRank.pres,
                                  n_common_genes = length(common))
  }
}

eff  <- do.call(rbind, res_effect)
pres <- do.call(rbind, res_pres)
# PATCH: keep earlier results for cohorts not re-run in this call
f_eff  <- file.path(out_dir, "01_ME11_effect_new_cohorts.csv")
f_pres <- file.path(out_dir, "01_module_preservation_new_cohorts.csv")
gse_of <- function(tag) sub("_.*", "", tag)
if (file.exists(f_eff)) {
  old <- read.csv(f_eff); old <- old[!gse_of(old$cohort) %in% cohorts, , drop = FALSE]
  if (nrow(old) && !is.null(eff)) old <- old[, colnames(eff)]
  eff <- rbind(old, eff)
}
if (file.exists(f_pres)) {
  old <- read.csv(f_pres); old <- old[!gse_of(old$cohort) %in% cohorts, , drop = FALSE]
  pres <- rbind(old, pres)
}
write.csv(eff,  file.path(out_dir, "01_ME11_effect_new_cohorts.csv"), row.names = FALSE)
write.csv(pres, file.path(out_dir, "01_module_preservation_new_cohorts.csv"), row.names = FALSE)

# (c) meta-analysis: original cohorts + new cohorts
# PATCH: the original cohorts are re-scored with the SAME ME11 score (mean z of the
# ME11 genes, z-scored within each dataset) instead of the paper's eigengene d values,
# so every cohort is measured the same way. Within a dataset, ComBat is a per-gene
# linear shift/scale, so z-scores are the same before and after ComBat.
orig <- do.call(rbind, lapply(c("GSE42133", "GSE111175"), function(s) {
  k  <- ref$pheno$dataset == s
  sc <- module_score(ref$Xgene[k, , drop = FALSE], ref$module_genes)
  cbind(cohort = s, cohens_d(sc, ref$pheno$diagnosis[k]))[, c("cohort", "d", "n_ASD", "n_TD", "var_d")]
}))
cat("Original cohorts, ME11 score (paper eigengene d: GSE42133 -0.79, GSE111175 +0.27):\n"); print(orig)
metas <- list()
for (def in c("broad", "strict")) {
  new  <- eff[eff$subset == paste0("all_", def), c("cohort", "d", "n_ASD", "n_TD", "var_d")]
  all_ <- rbind(orig, new)
  fe <- metafor::rma(yi = d, vi = var_d, data = all_, method = "FE")
  re <- metafor::rma(yi = d, vi = var_d, data = all_, method = "REML")
  # replication-only meta (new cohorts without the discovery cohorts)
  rep_ <- if (nrow(new) >= 2) metafor::rma(yi = d, vi = var_d, data = new, method = "REML") else NULL
  metas[[def]] <- data.frame(asd_definition = def, model = c("fixed (all)", "random (all)", "random (new only)"),
                     d = c(fe$b, re$b, if (is.null(rep_)) NA else rep_$b),
                     ci_lo = c(fe$ci.lb, re$ci.lb, if (is.null(rep_)) NA else rep_$ci.lb),
                     ci_hi = c(fe$ci.ub, re$ci.ub, if (is.null(rep_)) NA else rep_$ci.ub),
                     p = c(fe$pval, re$pval, if (is.null(rep_)) NA else rep_$pval),
                     Q = fe$QE, Q_p = fe$QEp, I2 = c(NA, re$I2, if (is.null(rep_)) NA else rep_$I2),
                     k = c(fe$k, re$k, nrow(new)))
  write.csv(all_, file.path(out_dir, paste0("01_meta_inputs_", def, ".csv")), row.names = FALSE)
  pdf(file.path(out_dir, paste0("01_ME11_forest_", def, ".pdf")), 7, 4)
  metafor::forest(re, slab = all_$cohort, header = paste("ME11 score, ASD vs TD (", def, ")")); dev.off()
}
meta <- do.call(rbind, metas)
write.csv(meta, file.path(out_dir, "01_ME11_meta_analysis.csv"), row.names = FALSE)
print(meta)

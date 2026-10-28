# =====================================================================
# 07_consensus_modules_replication.R   (new; not in the professor's set)
# Tests EVERY consensus module from 02 (built without ComBat) in the
# original and the new ASD cohorts, with the same score and labels as 01.
# Modules are compared on an equal footing, with FDR across modules,
# so that whichever module replicates is not cherry-picked.
#   Needs: results/02_consensus_labels.rds, geo_cache/ from 01
#   Outputs: results/07_consensus_by_cohort.csv, results/07_consensus_meta.csv
# =====================================================================
source("00_config.R")
need(c("GEOquery", "limma"), bioc = TRUE)
need(c("WGCNA", "metafor"))
ref <- load_reference()
cl  <- readRDS(file.path(out_dir, "02_consensus_labels.rds"))
sym <- ref$p2g$symbol[match(names(cl), ref$p2g$probe)]
ok  <- !is.na(sym) & sym != "" & !grepl("^LOC\\d+", sym) & cl != 0
sets <- lapply(split(sym[ok], cl[ok]), unique)
names(sets) <- paste0("cME", names(sets))
sets$ME11_original <- ref$module_genes
cat(length(sets), "gene sets\n")

# same helpers as 01 (copied so this script runs on its own)
src <- readLines("01_asd_replication.R")
grab <- function(fn) { i <- grep(paste0("^", fn, " <- function"), src); j <- i + which(src[i:length(src)] == "}")[1] - 1
  eval(parse(text = src[i:j]), envir = globalenv()) }
for (fn in c("guess_symbol", "guess_dx", "infer_sex")) grab(fn)

rows <- list()
add <- function(cohort, def, Xg, g) {
  for (m in names(sets)) {
    gg <- intersect(sets[[m]], colnames(Xg)); if (length(gg) < 10) next
    sc <- rowMeans(scale(Xg[, gg, drop = FALSE]), na.rm = TRUE)
    rows[[length(rows) + 1]] <<- cbind(cohort = cohort, asd_definition = def, module = m,
                                       n_genes = length(gg), cohens_d(sc, g))
  }
}
# original cohorts (z-scored within dataset, so ComBat does not matter)
for (s in c("GSE42133", "GSE111175")) for (def in c("broad", "strict")) {
  k <- ref$pheno$dataset == s; add(s, def, ref$Xgene[k, ], ref$pheno$diagnosis[k]) }

# new cohorts, from the files 01 already downloaded
for (gse in c("GSE18123", "GSE25507", "GSE26415")) {
  gl <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE, destdir = "geo_cache")
  for (es in gl) {
    tag <- paste0(gse, "_", Biobase::annotation(es))
    X <- t(Biobase::exprs(es)); if (max(X, na.rm = TRUE) > 100) X <- log2(pmax(X, 0) + 1)
    Xg <- collapse_max_mean(X, guess_symbol(Biobase::fData(es)))
    dx <- guess_dx(Biobase::pData(es))
    for (def in c("broad", "strict")) {
      g <- if (def == "broad") dx$dx else dx$dx_strict; k <- !is.na(g)
      if (length(table(g[k])) < 2 || min(table(g[k])) < 10) next
      add(tag, def, Xg[k, , drop = FALSE], g[k])
    }
    cat(tag, "done\n")
  }
}
byc <- do.call(rbind, rows)
write.csv(byc, file.path(out_dir, "07_consensus_by_cohort.csv"), row.names = FALSE)

disc <- c("GSE42133", "GSE111175")
meta <- do.call(rbind, lapply(split(byc, list(byc$module, byc$asd_definition), drop = TRUE), function(z) {
  new <- z[!z$cohort %in% disc, ]; old <- z[z$cohort %in% disc, ]
  fd <- metafor::rma(yi = d, vi = var_d, data = old, method = "FE")
  fn <- if (nrow(new) >= 2) metafor::rma(yi = d, vi = var_d, data = new, method = "REML") else NULL
  fa <- metafor::rma(yi = d, vi = var_d, data = z, method = "REML")
  data.frame(module = z$module[1], asd_definition = z$asd_definition[1], n_genes = z$n_genes[1],
             d_discovery = fd$b[1], p_discovery = fd$pval,
             d_replication = if (is.null(fn)) NA else fn$b[1], p_replication = if (is.null(fn)) NA else fn$pval,
             I2_replication = if (is.null(fn)) NA else fn$I2,
             n_new_same_direction = sum(sign(new$d) == sign(fd$b[1])), n_new = nrow(new),
             d_all = fa$b[1], p_all = fa$pval, I2_all = fa$I2)
}))
meta$FDR_discovery   <- ave(meta$p_discovery,   meta$asd_definition, FUN = function(p) p.adjust(p, "BH"))
meta$FDR_replication <- ave(meta$p_replication, meta$asd_definition, FUN = function(p) p.adjust(p, "BH"))
meta <- meta[order(meta$asd_definition, meta$p_replication), ]
write.csv(meta, file.path(out_dir, "07_consensus_meta.csv"), row.names = FALSE)
print(head(meta[meta$asd_definition == "strict", c("module","n_genes","d_discovery","FDR_discovery","d_replication","p_replication","FDR_replication","n_new_same_direction","n_new")], 12))

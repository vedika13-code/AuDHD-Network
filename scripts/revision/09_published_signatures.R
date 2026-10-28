# =====================================================================
# 09_published_signatures.R   (new)
# Do PUBLISHED ASD blood signatures replicate in independent cohorts?
#   Pramparo et al. 2015 (JAMA Psychiatry): DE genes, 4-module classifier genes,
#     and the 12 WGCNA modules (eSupplement 2)
#   Gazestani et al. 2019 (Nat Neurosci): DE genes (Table S2) and the
#     high-confidence DE-ASD network genes (Table S3)
# Signed sets: score = mean z(genes up in ASD) - mean z(genes down in ASD),
#   so a positive d means the published direction is reproduced.
# Unsigned sets (Pramparo modules): score = mean z of the module genes.
# The two original cohorts (GSE42133, GSE111175) overlap the data these
# signatures were derived from, so they are reported but NOT used as replication.
#   Needs: inputs/published_signatures.csv, geo_cache/ from 01
#   Outputs: results/09_signatures_by_cohort.csv, results/09_signatures_meta.csv
# =====================================================================
source("00_config.R")
need(c("GEOquery", "limma"), bioc = TRUE)
need(c("WGCNA", "metafor"))
ref <- load_reference()
sig <- read.csv("inputs/published_signatures.csv")
sets <- split(sig, sig$set)

src <- readLines("01_asd_replication.R")
grab <- function(fn) { i <- grep(paste0("^", fn, " <- function"), src); j <- i + which(src[i:length(src)] == "}")[1] - 1
  eval(parse(text = src[i:j]), envir = globalenv()) }
for (fn in c("guess_symbol", "guess_dx")) grab(fn)

score <- function(Xg, s) {
  Z <- scale(Xg)
  up <- intersect(s$gene[s$sign >= 0], colnames(Z)); dn <- intersect(s$gene[s$sign < 0], colnames(Z))
  if (all(s$sign == 0)) return(list(sc = rowMeans(Z[, up, drop = FALSE], na.rm = TRUE), n = length(up)))
  sc <- (if (length(up)) rowMeans(Z[, up, drop = FALSE], na.rm = TRUE) else 0) -
        (if (length(dn)) rowMeans(Z[, dn, drop = FALSE], na.rm = TRUE) else 0)
  list(sc = sc, n = length(up) + length(dn))
}
rows <- list()
add <- function(cohort, def, Xg, g) for (m in names(sets)) {
  r <- score(Xg, sets[[m]]); if (r$n < 10) next
  rows[[length(rows) + 1]] <<- cbind(cohort = cohort, asd_definition = def, signature = m,
                                     n_genes_total = nrow(sets[[m]]), n_genes_measured = r$n, cohens_d(r$sc, g))
}
# original cohorts: all probes, pre-ComBat, collapsed within dataset
Xfull <- readRDS(f_expr_precombat_full)
symf  <- ref$p2g$symbol[match(colnames(Xfull), ref$p2g$probe)]
for (s in c("GSE42133", "GSE111175")) {
  k <- ref$pheno$dataset == s; Xg <- collapse_max_mean(Xfull[ref$pheno$sample[k], , drop = FALSE], symf)
  for (def in c("broad", "strict")) add(s, def, Xg, ref$pheno$diagnosis[k])
}
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
write.csv(byc, file.path(out_dir, "09_signatures_by_cohort.csv"), row.names = FALSE)

orig <- c("GSE42133", "GSE111175")
meta <- do.call(rbind, lapply(split(byc, list(byc$signature, byc$asd_definition), drop = TRUE), function(z) {
  new <- z[!z$cohort %in% orig, ]; old <- z[z$cohort %in% orig, ]
  fo <- metafor::rma(yi = d, vi = var_d, data = old, method = "FE")
  fn <- metafor::rma(yi = d, vi = var_d, data = new, method = "REML")
  signed <- any(sets[[z$signature[1]]]$sign != 0)
  ref_sign <- if (signed) 1 else sign(fo$b[1])      # expected direction
  data.frame(signature = z$signature[1], asd_definition = z$asd_definition[1], signed = signed,
             d_original_cohorts = fo$b[1], p_original = fo$pval,
             d_replication = fn$b[1], ci_lo = fn$ci.lb, ci_hi = fn$ci.ub, p_replication = fn$pval, I2 = fn$I2,
             n_new_expected_direction = sum(sign(new$d) == ref_sign), n_new = nrow(new))
}))
meta$FDR_replication <- ave(meta$p_replication, meta$asd_definition, FUN = function(p) p.adjust(p, "BH"))
meta <- meta[order(meta$asd_definition, meta$p_replication), ]
write.csv(meta, file.path(out_dir, "09_signatures_meta.csv"), row.names = FALSE)
print(meta[meta$asd_definition == "strict", c("signature","signed","d_original_cohorts","d_replication","p_replication","FDR_replication","n_new_expected_direction","n_new")])

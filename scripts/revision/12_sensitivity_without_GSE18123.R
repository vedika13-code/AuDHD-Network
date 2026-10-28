# =====================================================================
# 12_sensitivity_without_GSE18123.R   (new)
# GSE18123 goes against the discovery direction for almost every signature.
# Re-estimate replication in the two other new cohorts (GSE25507, GSE26415) only.
#   Needs: results/09_signatures_by_cohort.csv, results/07_consensus_by_cohort.csv
#   Output: results/12_replication_without_GSE18123.csv
# =====================================================================
source("00_config.R"); need("metafor")
b  <- read.csv(file.path(out_dir, "09_signatures_by_cohort.csv"))
c7 <- read.csv(file.path(out_dir, "07_consensus_by_cohort.csv")); names(c7)[names(c7) == "module"] <- "signature"
cols <- c("cohort", "asd_definition", "signature", "d", "var_d")
x <- rbind(b[, cols], c7[c7$signature %in% c("ME11_original", "cME26"), cols])
x <- x[x$asd_definition == "strict" & !x$cohort %in% c("GSE42133", "GSE111175") & !grepl("GSE18123", x$cohort), ]
r <- do.call(rbind, lapply(split(x, x$signature), function(z) {
  f <- metafor::rma(yi = d, vi = var_d, data = z, method = "FE")
  data.frame(signature = z$signature[1], k = nrow(z), d = f$b[1], p = f$pval) }))
r$FDR <- p.adjust(r$p, "BH"); r <- r[order(r$p), ]
write.csv(r, file.path(out_dir, "12_replication_without_GSE18123.csv"), row.names = FALSE); print(r)

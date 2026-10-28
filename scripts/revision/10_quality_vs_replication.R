# =====================================================================
# 10_quality_vs_replication.R   (new)
# Across all consensus modules: does probe quality (mean intensity,
# share of predicted/pseudogene or unmapped probes) predict whether a
# module's ASD association replicates in the new cohorts?
#   Needs: results/02_consensus_labels.rds, results/07_consensus_meta.csv,
#          GPL10558 annotation from the main pipeline
#   Output: results/10_module_quality_vs_replication.csv
# =====================================================================
source("00_config.R")
suppressPackageStartupMessages(library(data.table))
gpl_file <- file.path(path.expand("~/Desktop/AuDHD_Network"), "data/raw/geo_cache/GPL10558.soft.gz")
if (!file.exists(gpl_file)) gpl_file <- Sys.getenv("GPL10558_SOFT")
ref <- load_reference()
cl  <- readRDS(file.path(out_dir, "02_consensus_labels.rds"))
pre <- readRDS(f_expr_precombat); mu <- colMeans(pre)
pct <- 100 * rank(mu) / length(mu); names(pct) <- colnames(pre)
g <- fread(gpl_file, skip = "ID\tSpecies", sep = "\t", fill = TRUE, quote = "")
a <- g[match(names(cl), g$ID)]
q <- data.frame(probe = names(cl), module = paste0("cME", cl), intensity_pct = pct[names(cl)],
                pred_pseudo = grepl("PREDICTED|pseudogene", a$Definition, ignore.case = TRUE),
                unmapped = is.na(a$Chromosome) | a$Chromosome == "",
                not_NM = !grepl("^NM_", a$RefSeq_ID))
q <- q[q$module != "cME0", ]
mq <- aggregate(cbind(intensity_pct, pred_pseudo, unmapped, not_NM) ~ module, q, mean)
m7 <- read.csv(file.path(out_dir, "07_consensus_meta.csv"))
m7 <- m7[m7$asd_definition == "strict" & m7$module != "ME11_original", ]
x <- merge(mq, m7, by = "module")
x$abs_disc <- abs(x$d_discovery)
x$retained <- x$d_replication * sign(x$d_discovery)        # replication effect in the discovery direction
x$shrink   <- x$retained / x$abs_disc
write.csv(x, file.path(out_dir, "10_module_quality_vs_replication.csv"), row.names = FALSE)
cat("Modules:", nrow(x), "\n")
for (v in c("intensity_pct", "pred_pseudo", "unmapped")) for (y in c("I2_replication", "retained", "abs_disc")) {
  ct <- cor.test(x[[v]], x[[y]], method = "spearman", exact = FALSE)
  cat(sprintf("%-14s vs %-15s rho = %5.2f  p = %.3g\n", v, y, ct$estimate, ct$p.value))
}

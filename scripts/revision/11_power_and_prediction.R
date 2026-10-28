# =====================================================================
# 11_power_and_prediction.R   (new)
# (a) Sample size per group needed to detect the discovery effect of ME11
#     and smaller effects (80% power, two-sided alpha 0.05).
# (b) Between-cohort heterogeneity (tau) of ME11 and the 95% prediction
#     interval for the effect in a NEW cohort.
#   Needs: results/01_meta_inputs_strict.csv, results/01_meta_inputs_broad.csv
#   Output: results/11_power_and_prediction.csv
# =====================================================================
source("00_config.R"); need("metafor")
n_needed <- function(d) ceiling(power.t.test(delta = d, sd = 1, power = 0.8, sig.level = 0.05)$n)
pw <- data.frame(effect_d = c(0.79, 0.55, 0.44, 0.38, 0.30, 0.20),
                 source = c("ME11 GSE42133 eigengene (paper)", "ME11 GSE42133 score", "ME11 discovery meta (score)",
                            "ME11 within-dataset meta (paper)", "", ""))
pw$n_per_group_80pct <- sapply(pw$effect_d, n_needed)
print(pw)
out <- list()
for (def in c("strict", "broad")) {
  x <- read.csv(file.path(out_dir, paste0("01_meta_inputs_", def, ".csv")))
  f <- metafor::rma(yi = d, vi = var_d, data = x, method = "REML")
  p <- predict(f)
  out[[def]] <- data.frame(asd_definition = def, k = f$k, d_pooled = f$b[1], tau = sqrt(f$tau2), I2 = f$I2,
                           pi_lo = p$pi.lb, pi_hi = p$pi.ub)
}
pr <- do.call(rbind, out); print(pr)
write.csv(merge(pw, pr, all = TRUE), file.path(out_dir, "11_power_and_prediction.csv"), row.names = FALSE)

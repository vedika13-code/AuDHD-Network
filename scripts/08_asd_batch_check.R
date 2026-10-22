# ============================================================
# Script: 08_asd_batch_check.R
# Purpose: Check whether the ME11-ASD association holds within
#          each source dataset (GSE111175, GSE42133) separately.
#
# Why: diagnosis is unbalanced across the two datasets
#      (GSE111175 ~28 ASD / 70 TD; GSE42133 ~91 ASD / 56 TD) and
#      ComBat was run with diagnosis in the model. With an
#      unbalanced design this can overstate group differences
#      (Nygaard et al., Biostatistics 2016). Within-dataset tests
#      are not affected, because ComBat applies the same
#      per-gene shift and scale to every sample in a dataset.
#
# Runs in about a minute. Output: results/review/batch_check.txt
# ============================================================

suppressPackageStartupMessages(library(WGCNA))
dir.create("results/review", recursive = TRUE, showWarnings = FALSE)
out <- "results/review/batch_check.txt"; cat("", file = out)
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = out, append = TRUE) }

asd      <- readRDS("data/processed/ASD_expression_clean.rds")
asd_meta <- readRDS("data/processed/ASD_metadata_clean.rds")
net_asd  <- readRDS("data/processed/ASD_network_final.rds")
asd_meta <- asd_meta[match(colnames(asd), rownames(asd_meta)), ]

ME11 <- moduleEigengenes(t(asd), net_asd$colors)$eigengenes[, "ME11"]
dx   <- factor(asd_meta$diagnosis, levels = c("TD", "ASD"))
bat  <- factor(asd_meta$batch)

say("Samples by dataset and diagnosis:")
capture.output(print(table(bat, dx)), file = out, append = TRUE); print(table(bat, dx))
say("Sex coding by dataset (GSE111175 uses 'M', GSE42133 uses 'male'):")
capture.output(print(table(bat, asd_meta$gender)), file = out, append = TRUE); print(table(bat, asd_meta$gender))

cohen_d <- function(x, g) {
  a <- x[g == "ASD"]; b <- x[g == "TD"]
  sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  d <- (mean(a) - mean(b)) / sp
  se <- sqrt(1 / length(a) + 1 / length(b) + d^2 / (2 * (length(a) + length(b))))
  c(d = d, se = se, lo = d - 1.96 * se, hi = d + 1.96 * se)
}

say("\n== ME11 eigengene, ASD vs TD ==")
rows <- list()
for (b in c("all", levels(bat))) {
  k <- if (b == "all") rep(TRUE, length(dx)) else bat == b
  f <- lm(ME11[k] ~ dx[k]); co <- summary(f)$coefficients[2, ]
  r <- cor(ME11[k], as.numeric(dx[k] == "ASD"))
  cd <- cohen_d(ME11[k], dx[k])
  rows[[b]] <- data.frame(subset = b, n = sum(k), n_ASD = sum(dx[k] == "ASD"),
                          r = r, d = cd["d"], d_lo = cd["lo"], d_hi = cd["hi"], p = co[4])
}
f_adj <- lm(ME11 ~ dx + bat); co <- summary(f_adj)$coefficients["dxASD", ]
tab <- do.call(rbind, rows); rownames(tab) <- NULL
print(tab); capture.output(print(tab), file = out, append = TRUE)
say("Pooled, dataset as covariate: beta = ", signif(co[1], 3), ", p = ", signif(co[4], 3))

# Fixed-effect meta-analysis of the two within-dataset effects (no reliance on ComBat)
w  <- tab[tab$subset != "all", ]
se <- (w$d_hi - w$d) / 1.96; wt <- 1 / se^2
d_meta <- sum(wt * w$d) / sum(wt); se_meta <- sqrt(1 / sum(wt))
p_meta <- 2 * pnorm(-abs(d_meta / se_meta))
Q <- sum(wt * (w$d - d_meta)^2); p_Q <- pchisq(Q, df = nrow(w) - 1, lower.tail = FALSE)
say("Meta-analysis of within-dataset effects: d = ", signif(d_meta, 3), " (95% CI ",
    signif(d_meta - 1.96 * se_meta, 3), " to ", signif(d_meta + 1.96 * se_meta, 3),
    "), p = ", signif(p_meta, 3))
say("Heterogeneity between datasets: Q = ", signif(Q, 3), ", p = ", signif(p_Q, 3))

# GSE111175 only: does age matter?
g111 <- read.csv("data/processed/GSE111175_metadata.csv", check.names = FALSE)
age_col <- grep("^age", names(g111), value = TRUE)[1]
age <- suppressWarnings(as.numeric(g111[[age_col]][match(rownames(asd_meta), g111$geo_accession)]))
k <- !is.na(age)
say("\nGSE111175 samples with age: ", sum(k))
for (fm in c("ME11 ~ dx", "ME11 ~ dx + age")) {
  dd <- data.frame(ME11 = ME11[k], dx = dx[k], age = age[k])
  co <- summary(lm(as.formula(fm), dd))$coefficients["dxASD", ]
  say(fm, ": beta = ", signif(co[1], 3), ", p = ", signif(co[4], 3))
}
say("Age by diagnosis (months): ",
    paste(names(tapply(age[k], dx[k], mean)), round(tapply(age[k], dx[k], mean), 1), collapse = "; "))

say("\nDone. Send results/review/batch_check.txt back.")

# ============================================================
# Script: 07_reviewer_analyses.R
# Purpose: Additional analyses requested in review.
#   A. ADHD limma accounting for twin pairing (+ inferred sex,
#      leukocyte marker scores as covariates)
#   B. ME11 module score: mixed model, paired test, Hedges' g,
#      TOST equivalence
#   C. Positive controls (score detects ME11 in ASD; spike-in
#      recovery in ADHD; sex-gene sanity check)
#   D. Leukocyte composition check (ASD and ADHD)
#   E. ASD module-trait model adjusted for sex / age
#   F. WGCNA modulePreservation (ASD reference -> ADHD test)
#   G. GO / Reactome re-run with the network gene universe
#
# Input:  data/processed/ASD_expression_clean.rds
#         data/processed/ASD_metadata_clean.rds
#         data/processed/ASD_network_final.rds
#         data/processed/ADHD_expression_clean.rds
#         data/processed/ADHD_metadata_clean.rds
#         data/processed/GSE159104_metadata.csv
#         data/processed/GSE111175_metadata.csv
#         data/processed/ME11_all_genes_annotated.csv
# Output: results/review/*.csv, results/review/summary.txt
#
# Run from the project root:  Rscript scripts/07_reviewer_analyses.R
# Section F (modulePreservation) takes the longest; set
# RUN_PRESERVATION <- FALSE to skip it on a first pass.
# ============================================================

suppressPackageStartupMessages({
  library(limma); library(WGCNA); library(lme4); library(lmerTest)
  library(clusterProfiler); library(ReactomePA); library(org.Hs.eg.db)
})
allowWGCNAThreads()
set.seed(2026)
RUN_PRESERVATION <- TRUE
dir.create("results/review", recursive = TRUE, showWarnings = FALSE)
log_file <- "results/review/summary.txt"
cat("", file = log_file)
say <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------
asd      <- readRDS("data/processed/ASD_expression_clean.rds")    # probes x samples
asd_meta <- readRDS("data/processed/ASD_metadata_clean.rds")
net_asd  <- readRDS("data/processed/ASD_network_final.rds")
adhd     <- readRDS("data/processed/ADHD_expression_clean.rds")   # genes x individuals
adhd_meta<- readRDS("data/processed/ADHD_metadata_clean.rds")
ME11_annot <- read.csv("data/processed/ME11_all_genes_annotated.csv")
ME11_sym <- unique(ME11_annot$Symbol)
ME11_sym <- ME11_sym[!grepl("^LOC", ME11_sym) & ME11_sym != "" & !is.na(ME11_sym)]

say("ADHD matrix as saved: ", nrow(adhd), " rows x ", ncol(adhd), " samples")

# The saved ADHD matrix (02_preprocessing.R) is indexed by UCSC transcript IDs
# (11,964 rows). The original analysis mapped these to 6,999 gene symbols, but the
# exact collapse step was never saved. So we rebuild it several ways, rerun the
# ORIGINAL limma model (~ condition, unpaired) for each, and keep the version whose
# t-statistics best reproduce results/adhd_de/ADHD_full_de_results.csv.
if (mean(grepl("^uc", rownames(adhd))) > 0.5) {
  count_file <- list.files("data/raw/GSE159104", pattern = "DGE_raw", full.names = TRUE)[1]
  stopifnot(!is.na(count_file))
  ann <- read.table(gzfile(count_file), header = TRUE, row.names = 1,
                    sep = "\t", check.names = FALSE)[, 1:7]          # same call as 02_preprocessing.R
  sym_col <- grep("gene.?name", names(ann), ignore.case = TRUE, value = TRUE)[1]
  say("Annotation column used for symbols: ", sym_col)
  sym <- as.character(ann[[sym_col]][match(rownames(adhd), rownames(ann))])
  ok  <- !is.na(sym) & sym != ""
  X <- adhd[ok, ]; sym <- sym[ok]

  collapse <- list(
    mean    = function(M) limma::avereps(M, ID = sym),
    first   = function(M) { k <- !duplicated(sym); r <- M[k, ]; rownames(r) <- sym[k]; r },
    maxmean = function(M) { o <- order(-rowMeans(M, na.rm = TRUE)); k <- o[!duplicated(sym[o])]
                            r <- M[k, ]; rownames(r) <- sym[k]; r }
  )
  transform <- list(
    as_saved = function(M) M,
    log2p1   = function(M) log2(M + 1),
    rowz     = function(M) t(scale(t(M)))
  )
  saved <- read.csv("results/adhd_de/ADHD_full_de_results.csv", row.names = 1)
  cond0 <- factor(adhd_meta$ADHD_Status[match(colnames(adhd), rownames(adhd_meta))], levels = c("CTRL", "ADHD"))
  d0 <- model.matrix(~ cond0)
  cmp <- list(); mats <- list()
  for (cn in names(collapse)) for (tn in names(transform)) {
    M <- transform[[tn]](collapse[[cn]](X))
    M <- M[rowSums(is.na(M)) == 0, , drop = FALSE]
    tt <- topTable(eBayes(lmFit(M, d0)), coef = 2, number = Inf, sort.by = "none")
    g <- intersect(rownames(tt), rownames(saved))
    key <- paste(cn, tn, sep = "+"); mats[[key]] <- M
    cmp[[key]] <- data.frame(version = key, n_genes = nrow(M), n_shared = length(g),
                             r_t = cor(tt[g, "t"], saved[g, "t"]),
                             r_logFC = cor(tt[g, "logFC"], saved[g, "logFC"]))
  }
  cmp <- do.call(rbind, cmp); cmp <- cmp[order(-cmp$r_t), ]
  print(cmp); capture.output(print(cmp), file = log_file, append = TRUE)
  write.csv(cmp, "results/review/ADHD_collapse_reconstruction.csv", row.names = FALSE)
  best <- cmp$version[1]
  adhd <- mats[[best]]
  say("Using '", best, "': ", nrow(adhd), " genes; reproduces original t-statistics with r = ",
      round(cmp$r_t[1], 4))
  if (cmp$r_t[1] < 0.99) say("WARNING: no version reproduced the original DE closely (r < 0.99). ",
                             "Send results/review/ADHD_collapse_reconstruction.csv back before using these results.")
}

# ------------------------------------------------------------
# Twin pair IDs from the GEO series metadata
# ------------------------------------------------------------
geo <- read.csv("data/processed/GSE159104_metadata.csv", check.names = FALSE)
geo$individual <- sub("_rep[0-9]+$", "", geo$title)
geo$pair <- trimws(sub("twin pair id:", "", geo[["twin pair id:ch1"]]))
pair_map <- unique(geo[, c("individual", "pair")])
pair_map <- pair_map[!duplicated(pair_map$individual), ]

ids <- colnames(adhd)
adhd_meta <- adhd_meta[match(ids, rownames(adhd_meta)), ]
stopifnot(all(rownames(adhd_meta) == ids))
adhd_meta$pair <- pair_map$pair[match(ids, pair_map$individual)]
adhd_meta$block <- ifelse(is.na(adhd_meta$pair) | adhd_meta$pair %in% c("NA", ""),
                          ids, paste0("pair_", adhd_meta$pair))
adhd_meta$condition <- factor(adhd_meta$ADHD_Status, levels = c("CTRL", "ADHD"))
say("Twin pairs: ", sum(table(adhd_meta$block) == 2), "; unrelated individuals: ",
    sum(table(adhd_meta$block) == 1))

# ------------------------------------------------------------
# Inferred sex (no sex/age in GEO metadata for GSE159104)
# ------------------------------------------------------------
y_genes <- intersect(c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY"), rownames(adhd))
say("Y-linked genes available for sex inference: ", paste(y_genes, collapse = ", "))
y_score <- colMeans(t(scale(t(adhd[y_genes, , drop = FALSE]))))
km <- kmeans(y_score, centers = 2)
male_cluster <- which.max(tapply(y_score, km$cluster, mean))
adhd_meta$sex_inferred <- factor(ifelse(km$cluster == male_cluster, "M", "F"))
say("Inferred sex: ", paste(names(table(adhd_meta$sex_inferred)), table(adhd_meta$sex_inferred), collapse = "; "))
say("Sex x diagnosis: "); capture.output(print(table(adhd_meta$sex_inferred, adhd_meta$condition)), file = log_file, append = TRUE)
write.csv(data.frame(sample = ids, y_score, sex_inferred = adhd_meta$sex_inferred),
          "results/review/ADHD_inferred_sex.csv", row.names = FALSE)

# ------------------------------------------------------------
# D. Leukocyte marker scores (marker-based proxy for composition)
# ------------------------------------------------------------
markers <- list(
  neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "MNDA", "S100A12"),
  monocyte   = c("CD14", "LYZ", "FCN1", "CSF1R", "VCAN", "S100A9"),
  Bcell      = c("MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "CD22"),
  Tcell      = c("CD3E", "CD3D", "CD3G", "CD2", "CD6", "LCK"),
  NKcell     = c("GNLY", "NKG7", "KLRD1", "PRF1", "GZMB", "KLRF1"),
  eosinophil = c("CLC", "IL5RA", "CCR3", "SIGLEC8")
)
cell_score <- function(mat, symbols_of_rows) {
  out <- sapply(names(markers), function(ct) {
    g <- which(symbols_of_rows %in% markers[[ct]])
    if (length(g) < 2) return(rep(NA_real_, ncol(mat)))
    colMeans(t(scale(t(mat[g, , drop = FALSE]))), na.rm = TRUE)
  })
  rownames(out) <- colnames(mat); out
}
cs_adhd <- cell_score(adhd, rownames(adhd))
say("\n== D. Leukocyte marker scores, ADHD vs CTRL (mixed model, twin pair random effect) ==")
cell_res_adhd <- do.call(rbind, lapply(colnames(cs_adhd), function(ct) {
  if (all(is.na(cs_adhd[, ct]))) return(data.frame(cell = ct, estimate = NA, p = NA, n_markers = 0))
  d <- data.frame(y = cs_adhd[, ct], condition = adhd_meta$condition, block = adhd_meta$block)
  m <- lmer(y ~ condition + (1 | block), data = d)
  co <- summary(m)$coefficients["conditionADHD", ]
  data.frame(cell = ct, estimate = co["Estimate"], p = co["Pr(>|t|)"],
             n_markers = sum(rownames(adhd) %in% markers[[ct]]))
}))
print(cell_res_adhd); write.csv(cell_res_adhd, "results/review/ADHD_cell_scores_test.csv", row.names = FALSE)
capture.output(print(cell_res_adhd), file = log_file, append = TRUE)
cell_cov_adhd <- cs_adhd[, colSums(is.na(cs_adhd)) == 0, drop = FALSE]

# ASD side: probe -> symbol
annot_full <- tryCatch({
  library(GEOquery); Table(getGEO("GPL10558"))[, c("ID", "Symbol", "Entrez_Gene_ID")]
}, error = function(e) NULL)
stopifnot(!is.null(annot_full))
asd_sym <- annot_full$Symbol[match(rownames(asd), annot_full$ID)]
cs_asd <- cell_score(asd, asd_sym)
asd_meta <- asd_meta[match(colnames(asd), rownames(asd_meta)), ]
diag_bin <- as.numeric(asd_meta$diagnosis == "ASD")

MEs <- moduleEigengenes(t(asd), net_asd$colors)$eigengenes
ME11 <- MEs[, "ME11"]
say("\n== D. ASD: cell scores vs diagnosis, and vs ME11 eigengene ==")
cell_res_asd <- do.call(rbind, lapply(colnames(cs_asd), function(ct) {
  x <- cs_asd[, ct]; if (all(is.na(x))) return(NULL)
  data.frame(cell = ct,
             p_vs_diagnosis = t.test(x ~ asd_meta$diagnosis)$p.value,
             r_with_ME11 = cor(x, ME11), p_r = cor.test(x, ME11)$p.value)
}))
print(cell_res_asd); capture.output(print(cell_res_asd), file = log_file, append = TRUE)
write.csv(cell_res_asd, "results/review/ASD_cell_scores_test.csv", row.names = FALSE)

# ------------------------------------------------------------
# E. ASD module-trait adjusted for sex, (age), cell scores
# ------------------------------------------------------------
say("\n== E. ME11 ~ diagnosis, adjusted ==")
g111 <- read.csv("data/processed/GSE111175_metadata.csv", check.names = FALSE)
age_col <- grep("^age", names(g111), value = TRUE)[1]
asd_meta$age_months <- suppressWarnings(as.numeric(g111[[age_col]][match(rownames(asd_meta), g111$geo_accession)]))
say("Age available for ", sum(!is.na(asd_meta$age_months)), " of ", nrow(asd_meta), " ASD-cohort samples")
dA <- data.frame(ME11, dx = factor(asd_meta$diagnosis, levels = c("TD", "ASD")),
                 sex = factor(asd_meta$gender), cs_asd[, colSums(is.na(cs_asd)) == 0, drop = FALSE],
                 age = asd_meta$age_months)
fits <- list(
  unadjusted   = lm(ME11 ~ dx, dA),
  sex          = lm(ME11 ~ dx + sex, dA),
  sex_cells    = lm(as.formula(paste("ME11 ~ dx + sex +", paste(colnames(cs_asd)[colSums(is.na(cs_asd)) == 0], collapse = "+"))), dA),
  sex_age_GSE111175_only = lm(ME11 ~ dx + sex + age, dA[!is.na(dA$age), ])
)
adj_tab <- do.call(rbind, lapply(names(fits), function(n) {
  co <- summary(fits[[n]])$coefficients["dxASD", ]
  data.frame(model = n, n = nobs(fits[[n]]), beta = co[1], se = co[2], p = co[4])
}))
print(adj_tab); capture.output(print(adj_tab), file = log_file, append = TRUE)
write.csv(adj_tab, "results/review/ASD_ME11_adjusted.csv", row.names = FALSE)

# ------------------------------------------------------------
# A. ADHD limma with twin pairing (duplicateCorrelation)
# ------------------------------------------------------------
say("\n== A. limma, twin pair as blocking factor ==")
run_limma <- function(design, label) {
  dc  <- duplicateCorrelation(adhd, design, block = adhd_meta$block)
  fit <- eBayes(lmFit(adhd, design, block = adhd_meta$block, correlation = dc$consensus))
  tt  <- topTable(fit, coef = "conditionADHD", number = Inf, sort.by = "P")
  say(label, ": consensus within-pair correlation = ", round(dc$consensus, 3),
      "; FDR<0.05: ", sum(tt$adj.P.Val < 0.05), "; FDR<0.10: ", sum(tt$adj.P.Val < 0.10))
  tt
}
des_pair <- model.matrix(~ condition, adhd_meta)
tt_pair  <- run_limma(des_pair, "Paired (block)")
des_cov  <- model.matrix(~ condition + sex_inferred + cell_cov_adhd, adhd_meta)
colnames(des_cov) <- make.names(colnames(des_cov))
tt_cov   <- run_limma(des_cov, "Paired + sex + cell scores")
write.csv(tt_pair, "results/review/ADHD_limma_paired.csv")
write.csv(tt_cov,  "results/review/ADHD_limma_paired_covariates.csv")

# Fisher test for ME11 among nominal DE genes (FDR gives 0 genes, so use p<0.05)
ME11_in <- intersect(ME11_sym, rownames(adhd))
fisher_nominal <- function(tt) {
  de <- rownames(tt)[tt$P.Value < 0.05]; N <- nrow(tt); k <- length(intersect(ME11_in, de))
  ft <- fisher.test(matrix(c(k, length(ME11_in) - k, length(de) - k, N - length(ME11_in) - length(de) + k), 2),
                    alternative = "greater")
  c(overlap = k, n_ME11 = length(ME11_in), n_de = length(de), OR = unname(ft$estimate), p = ft$p.value)
}
say("Fisher (nominal p<0.05), paired: ", paste(names(fisher_nominal(tt_pair)), signif(fisher_nominal(tt_pair), 3), collapse = " "))
say("Fisher (nominal p<0.05), paired+cov: ", paste(names(fisher_nominal(tt_cov)), signif(fisher_nominal(tt_cov), 3), collapse = " "))

# Competitive gene-set test that respects the block correlation
cam <- camera(adhd, index = which(rownames(adhd) %in% ME11_in), design = des_pair,
              contrast = "conditionADHD")
say("camera (ME11 set, ADHD vs CTRL): ", paste(names(cam), unlist(cam), collapse = " "))

# ------------------------------------------------------------
# B. Module score: mixed model, paired test, Hedges' g, TOST
# ------------------------------------------------------------
say("\n== B. ME11 module score in ADHD ==")
mod_score <- function(mat, genes) colMeans(t(scale(t(mat[genes, , drop = FALSE]))), na.rm = TRUE)
adhd_meta$score <- mod_score(adhd, ME11_in)
sdT <- sd(adhd_meta$score)

m0 <- lmer(score ~ condition + (1 | block), data = adhd_meta)
cc <- as.data.frame(cell_cov_adhd); names(cc) <- paste0("cell_", names(cc))
adhd_meta <- cbind(adhd_meta, cc)
m1 <- lmer(as.formula(paste("score ~ condition + sex_inferred +", paste(names(cc), collapse = " + "), "+ (1 | block)")),
           data = adhd_meta)
for (nm in c("m0", "m1")) {
  co <- summary(get(nm))$coefficients["conditionADHD", ]
  ci <- co["Estimate"] + c(-1, 1) * qt(0.975, co["df"]) * co["Std. Error"]
  say(nm, ": beta = ", signif(co["Estimate"], 3), " (", signif(co["Estimate"] / sdT, 3), " SD), 95% CI ",
      signif(ci[1] / sdT, 3), " to ", signif(ci[2] / sdT, 3), " SD, p = ", signif(co["Pr(>|t|)"], 3))
  # TOST at +/- 0.5 SD
  M <- 0.5 * sdT; se <- co["Std. Error"]; df <- co["df"]; b <- co["Estimate"]
  p_tost <- max(pt((b + M) / se, df, lower.tail = FALSE), pt((b - M) / se, df))
  say(nm, ": TOST (+/-0.5 SD) p = ", signif(p_tost, 3))
}
tw <- adhd_meta[grepl("^pair_", adhd_meta$block), ]
pw <- reshape(tw[, c("block", "condition", "score")], idvar = "block", timevar = "condition", direction = "wide")
pt_res <- t.test(pw$score.ADHD, pw$score.CTRL, paired = TRUE)
say("Paired t-test (", nrow(pw), " discordant pairs): p = ", signif(pt_res$p.value, 3))
a <- adhd_meta$score[adhd_meta$condition == "ADHD"]; c0 <- adhd_meta$score[adhd_meta$condition == "CTRL"]
n1 <- length(a); n2 <- length(c0); sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(c0)) / (n1 + n2 - 2))
g <- (mean(a) - mean(c0)) / sp * (1 - 3 / (4 * (n1 + n2) - 9)); seg <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)))
say("Hedges' g = ", signif(g, 3), " (95% CI ", signif(g - 1.96 * seg, 3), " to ", signif(g + 1.96 * seg, 3), ")")

# ------------------------------------------------------------
# C. Positive controls
# ------------------------------------------------------------
say("\n== C. Positive controls ==")
# C1. Same scoring procedure detects ME11 in the ASD cohort (restricted to the
#     same genes that are available in the ADHD matrix, so the gene set is matched)
asd_probe_ME11 <- names(net_asd$colors)[net_asd$colors == 11]
asd_probe_matched <- asd_probe_ME11[asd_sym[match(asd_probe_ME11, rownames(asd))] %in% ME11_in]
sA <- mod_score(asd, asd_probe_matched)
xA <- sA[asd_meta$diagnosis == "ASD"]; yA <- sA[asd_meta$diagnosis == "TD"]
spA <- sqrt(((length(xA) - 1) * var(xA) + (length(yA) - 1) * var(yA)) / (length(xA) + length(yA) - 2))
say("C1. ASD cohort, matched-gene ME11 score (", length(asd_probe_matched), " probes): d = ",
    signif((mean(xA) - mean(yA)) / spA, 3), ", p = ", signif(t.test(xA, yA)$p.value, 3))

# C2. Spike-in recovery: add a shift of d = 0.54 (ASD effect size) to ME11 genes in
#     ADHD cases, recompute the score and the mixed model; report detection rate.
#     (Shift applied on the module-score scale, labels permuted each time so any
#      true signal is removed and only the spiked effect remains.)
spike_score <- function(d_target, nrep = 500) {
  mean(replicate(nrep, {
    perm <- sample(adhd_meta$condition)               # break any true signal
    s <- adhd_meta$score + d_target * sdT * (perm == "ADHD")
    d <- data.frame(s, perm, block = adhd_meta$block)
    summary(lmer(s ~ perm + (1 | block), d))$coefficients[2, "Pr(>|t|)"] < 0.05
  }))
}
for (dd in c(0.3, 0.54, 0.8)) say("C2. Spike-in d = ", dd, " on score scale: detection rate = ", spike_score(dd))

# C3. Sex-inference check: the Y-linked score should be clearly bimodal.
say("C3. Y-score gap between inferred sexes (SD units): ",
    signif(diff(tapply(y_score, adhd_meta$sex_inferred, mean)) / sd(y_score), 3))

# ------------------------------------------------------------
# F. Module preservation (ASD reference -> ADHD test)
# ------------------------------------------------------------
if (RUN_PRESERVATION) {
  say("\n== F. modulePreservation ==")
  # Collapse ASD probes to gene symbols (max mean expression probe per gene)
  keep <- !is.na(asd_sym) & asd_sym != ""
  cr <- collapseRows(asd[keep, ], rowGroup = asd_sym[keep], rowID = rownames(asd)[keep],
                     method = "MaxMean")
  asd_gene <- cr$datETcollapsed
  col_gene <- net_asd$colors[cr$group2row[, "selectedRowID"]]
  names(col_gene) <- cr$group2row[, "group"]
  common <- intersect(rownames(asd_gene), rownames(adhd))
  say("Genes shared between ASD (collapsed) and ADHD: ", length(common),
      "; ME11 genes among them: ", sum(col_gene[common] == 11))
  multiExpr  <- list(ASD = list(data = t(asd_gene[common, ])), ADHD = list(data = t(adhd[common, ])))
  multiColor <- list(ASD = col_gene[common])
  mp <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1, nPermutations = 200,
                           networkType = "signed", randomSeed = 2026, quickCor = 0,
                           maxModuleSize = 1000, verbose = 3)
  st <- mp$preservation$Z$ref.ASD$inColumnsAlsoPresentIn.ADHD
  ob <- mp$preservation$observed$ref.ASD$inColumnsAlsoPresentIn.ADHD
  pres <- data.frame(module = rownames(st), moduleSize = st$moduleSize,
                     Zsummary = st$Zsummary.pres, medianRank = ob$medianRank.pres)
  print(pres); capture.output(print(pres), file = log_file, append = TRUE)
  write.csv(pres, "results/review/modulePreservation_ASD_to_ADHD.csv", row.names = FALSE)
  say("ME11 Zsummary = ", signif(pres$Zsummary[pres$module == "11"], 3),
      " (Z<2 no evidence, 2-10 weak/moderate, >10 strong preservation)")
}

# ------------------------------------------------------------
# G. Enrichment with the network universe
# ------------------------------------------------------------
say("\n== G. Enrichment, universe = all genes in the 15,000-probe ASD network ==")
universe <- unique(as.character(annot_full$Entrez_Gene_ID[match(rownames(asd), annot_full$ID)]))
universe <- universe[!is.na(universe) & universe != ""]
ME11_entrez <- unique(as.character(ME11_annot$Entrez_Gene_ID[!grepl("^LOC", ME11_annot$Symbol)]))
ME11_entrez <- ME11_entrez[!is.na(ME11_entrez) & ME11_entrez != ""]
say("Universe size: ", length(universe), "; ME11 genes: ", length(ME11_entrez))
for (ont in c("BP", "MF", "CC")) {
  r <- enrichGO(ME11_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
                universe = universe, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  write.csv(r@result, paste0("results/review/ME11_GO_", ont, "_networkUniverse.csv"), row.names = FALSE)
  say("GO ", ont, ": FDR<0.05 terms = ", sum(r@result$p.adjust < 0.05))
}
rp <- enrichPathway(ME11_entrez, organism = "human", universe = universe,
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
write.csv(rp@result, "results/review/ME11_Reactome_networkUniverse.csv", row.names = FALSE)
say("Reactome: FDR<0.05 pathways = ", sum(rp@result$p.adjust < 0.05))
capture.output(print(head(rp@result[, c("Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust")], 10)),
               file = log_file, append = TRUE)

say("\nDone. Send results/review/summary.txt back for the manuscript update.")

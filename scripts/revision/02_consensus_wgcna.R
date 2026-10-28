# =====================================================================
# 02_consensus_wgcna.R
# Build modules WITHOUT ComBat: a consensus network across GSE111175 and
# GSE42133, then test each consensus module within each dataset and
# meta-analyse. Also maps consensus modules back to the original ME11.
# Fixes reviewer point 2 (unbalanced ComBat, grey-module signal).
# Outputs: results/02_*.csv
# =====================================================================
source("00_config.R")
need(c("WGCNA", "metafor"))
allowWGCNAThreads()

X     <- readRDS(f_expr_precombat)
pheno <- readRDS(f_pheno)
pheno <- pheno[match(rownames(X), pheno$sample), ]
labs  <- readRDS(f_labels)[colnames(X)]
sets  <- unique(pheno$dataset)

multiExpr <- lapply(sets, function(s) list(data = X[pheno$dataset == s, , drop = FALSE]))
names(multiExpr) <- sets
gsg <- goodSamplesGenesMS(multiExpr, verbose = 0)
if (!gsg$allOK) for (s in seq_along(multiExpr))
  multiExpr[[s]]$data <- multiExpr[[s]]$data[gsg$goodSamples[[s]], gsg$goodGenes]

# Soft-threshold check per dataset (report this in the Methods)
sft <- do.call(rbind, lapply(sets, function(s) {
  r <- pickSoftThreshold(multiExpr[[s]]$data, powerVector = 1:20,
                         networkType = "signed", verbose = 0)$fitIndices
  cbind(dataset = s, r[, c("Power", "SFT.R.sq", "slope", "mean.k.")])
}))
write.csv(sft, file.path(out_dir, "02_soft_threshold_by_dataset.csv"), row.names = FALSE)

net <- blockwiseConsensusModules(
  multiExpr, power = soft_power, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE, numericLabels = TRUE, maxBlockSize = max_block,   # PATCH: was 20000 (one block, ~6 GB RAM)
  saveTOMs = FALSE, verbose = 1)
cons_labels <- setNames(net$colors, colnames(multiExpr[[1]]$data))
saveRDS(cons_labels, file.path(out_dir, "02_consensus_labels.rds"))
cat("Consensus modules:", length(setdiff(unique(cons_labels), 0)),
    " grey fraction:", round(mean(cons_labels == 0), 3), "\n")

# Within-dataset association of each consensus eigengene with diagnosis
rows <- list()
for (s in seq_along(sets)) {
  ME <- net$multiMEs[[s]]$data
  smp <- rownames(multiExpr[[s]]$data)
  dx  <- pheno$diagnosis[match(smp, pheno$sample)]
  for (m in colnames(ME)) rows[[length(rows) + 1]] <-
    cbind(dataset = sets[s], module = m, cohens_d(ME[, m], dx))
}
w <- do.call(rbind, rows)
meta <- do.call(rbind, lapply(split(w, w$module), function(z) {
  f <- metafor::rma(yi = z$d, vi = z$var_d, method = "FE")
  data.frame(module = z$module[1], d_meta = as.numeric(f$b), ci_lo = f$ci.lb, ci_hi = f$ci.ub,
             p = f$pval, Q = f$QE, Q_p = f$QEp,
             d_GSE111175 = z$d[z$dataset == "GSE111175"], d_GSE42133 = z$d[z$dataset == "GSE42133"])
}))
meta$FDR <- p.adjust(meta$p, "BH")
meta <- meta[order(meta$p), ]
write.csv(w,    file.path(out_dir, "02_consensus_module_by_dataset.csv"), row.names = FALSE)
write.csv(meta, file.path(out_dir, "02_consensus_module_meta.csv"), row.names = FALSE)

# Which consensus module(s) correspond to the original ME11?
ot <- overlapTable(cons_labels, labs[names(cons_labels)])
col <- as.character(target_module)
match_tab <- data.frame(consensus_module = rownames(ot$countTable),
                        overlap_with_ME11 = ot$countTable[, col],
                        p_overlap = ot$pTable[, col])
match_tab <- match_tab[order(match_tab$p_overlap), ]
match_tab$ME11_genes_total <- sum(labs == target_module)
write.csv(match_tab, file.path(out_dir, "02_consensus_vs_ME11_overlap.csv"), row.names = FALSE)
print(head(meta, 10)); print(head(match_tab, 5))

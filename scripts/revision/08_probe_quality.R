# 08_probe_quality.R  (new): probe quality of ME11 vs other modules (GPL10558 annotation),
# and whether the ME11-ASD association depends on which probes are used.
# Outputs: results/08_probe_quality_by_group.csv, results/08_ME11_by_probe_subset.csv
suppressMessages({library(WGCNA); library(data.table)})
suppressMessages(source("00_config.R")); ref <- load_reference()
g <- fread(file.path(path.expand("~/Desktop/AuDHD_Network"), "data/raw/geo_cache/GPL10558.soft.gz"), skip="ID\tSpecies", sep="\t", fill=TRUE, quote="")
g <- g[ID %in% colnames(ref$expr)]
a <- g[match(colnames(ref$expr), g$ID)]
labs <- ref$labels; is11 <- labs == 11
pre <- readRDS("inputs/asd_expr_precombat.rds")
mu  <- colMeans(pre)                         # mean log2 intensity before ComBat
sym <- a$Symbol
# primary probe = highest-mean probe of each symbol among ALL 47k probes
full <- readRDS("inputs/asd_expr_precombat_full.rds"); mfull <- colMeans(full)
p2g <- ref$p2g; sfull <- p2g$symbol[match(colnames(full), p2g$probe)]
prim <- tapply(seq_along(mfull), sfull, function(i) colnames(full)[i[which.max(mfull[i])]])
is_primary <- colnames(ref$expr) %in% unlist(prim)
nprobes <- table(sfull)[sym]; nprobes[is.na(nprobes)] <- 1
qual <- data.frame(
  group = ifelse(is11, "ME11", ifelse(labs == 0, "grey", "other modules")),
  mean_intensity_pct = 100 * rank(mu) / length(mu),
  secondary_probe = !is_primary & !is.na(sym) & sym != "",
  no_symbol_or_LOC = is.na(sym) | sym == "" | grepl("^LOC\\d+", sym),
  unmapped = is.na(a$Chromosome) | a$Chromosome == "",
  not_RefSeq_NM = !grepl("^NM_", a$RefSeq_ID),
  predicted_or_pseudo = grepl("PREDICTED|pseudogene", a$Definition, ignore.case = TRUE))
summ <- aggregate(. ~ group, qual, function(x) round(mean(x), 3))
cat("\n== Probe quality by group (fraction, or mean percentile) ==\n"); print(summ)
cat("\nMedian intensity percentile: "); print(tapply(qual$mean_intensity_pct, qual$group, median))
cat("Wilcoxon ME11 vs other modules intensity p =", signif(wilcox.test(qual$mean_intensity_pct[qual$group=="ME11"], qual$mean_intensity_pct[qual$group=="other modules"])$p.value,3), "\n")
for (v in c("secondary_probe","no_symbol_or_LOC","unmapped","not_RefSeq_NM","predicted_or_pseudo")) {
  t <- table(qual$group[qual$group!="grey"]=="ME11", qual[[v]][qual$group!="grey"])
  cat(v, ": Fisher p =", signif(fisher.test(t)$p.value,3), "\n") }
# Does the ASD association come from primary or secondary probes?
cat("\n== ME11 association by probe subset (eigengene d within dataset) ==\n")
res <- list()
for (sub in c("all ME11 probes","primary/unique probes","secondary probes","high intensity (top half)","low intensity (bottom half)")) {
  k <- switch(sub, "all ME11 probes"=is11, "primary/unique probes"=is11 & !qual$secondary_probe,
              "secondary probes"=is11 & qual$secondary_probe,
              "high intensity (top half)"=is11 & mu >= median(mu[is11]), "low intensity (bottom half)"=is11 & mu < median(mu[is11]))
  for (s in c("GSE42133","GSE111175")) {
    i <- ref$pheno$dataset == s; X <- ref$expr[i, k]
    me <- moduleEigengenes(X, rep(1, ncol(X)))$eigengenes[,1]
    if (cor(me, rowMeans(scale(X))) < 0) me <- -me
    cd <- cohens_d(me, ref$pheno$diagnosis[i])
    res[[length(res)+1]] <- data.frame(subset=sub, dataset=s, n_probes=sum(k), d=round(cd$d,2), ci_lo=round(cd$ci_lo,2), ci_hi=round(cd$ci_hi,2))
  } }
r <- do.call(rbind, res); print(r)
write.csv(summ, "results/08_probe_quality_by_group.csv", row.names=FALSE)
write.csv(r, "results/08_ME11_by_probe_subset.csv", row.names=FALSE)

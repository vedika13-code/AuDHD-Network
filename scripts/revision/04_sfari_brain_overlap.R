# =====================================================================
# 04_sfari_brain_overlap.R
# Anchor the blood modules to ASD genetics (SFARI) and ASD brain data
# (Gandal et al. 2022 cortex DE). Background = genes in the network.
# Fixes reviewer points 2 and 3 (biological relevance / take-home).
#
# Inputs you download yourself (see README):
#   inputs/SFARI-Gene_genes.csv      from gene.sfari.org (columns incl. gene-symbol, gene-score, syndromic)
#   inputs/gandal2022_cortex_DE.csv  two-three columns renamed to: gene, logFC, FDR
# Outputs: results/04_*.csv
# =====================================================================
source("00_config.R")
ref  <- load_reference()
bg   <- names(ref$gene_labels)                       # network genes (symbols)
mods <- ref$module_sets                                # PATCH: all probes of each module (ME11 = 187 genes)

fisher_row <- function(set, hits, bg, label) {
  set <- intersect(set, bg); hits <- intersect(hits, bg)
  a <- length(intersect(set, hits)); b <- length(set) - a
  c <- length(hits) - a;             d <- length(bg) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater")
  data.frame(test = label, module_size = length(set), overlap = a,
             expected = round(length(set) * length(hits) / length(bg), 2),
             OR = unname(ft$estimate), p = ft$p.value,
             genes = paste(sort(intersect(set, hits)), collapse = ";"))
}

## ---- SFARI ---------------------------------------------------------------
sf <- read.csv("inputs/SFARI-Gene_genes.csv", check.names = FALSE)
colnames(sf) <- tolower(gsub("[^A-Za-z]", "", colnames(sf)))      # genesymbol, genescore, syndromic
sf_all <- unique(sf$genesymbol)
sf_hc  <- unique(sf$genesymbol[sf$genescore %in% c(1, "1") | sf$syndromic %in% c(1, "1")])
sfari <- do.call(rbind, lapply(names(mods), function(m) rbind(
  cbind(module = m, fisher_row(mods[[m]], sf_all, bg, "SFARI all")),
  cbind(module = m, fisher_row(mods[[m]], sf_hc,  bg, "SFARI score 1 or syndromic")))))
sfari$FDR <- ave(sfari$p, sfari$test, FUN = function(p) p.adjust(p, "BH"))
write.csv(sfari[order(sfari$p), ], file.path(out_dir, "04_SFARI_enrichment_all_modules.csv"), row.names = FALSE)

## ---- ASD cortex (Gandal 2022) -------------------------------------------
if (file.exists("inputs/gandal2022_cortex_DE.csv")) {
  br <- read.csv("inputs/gandal2022_cortex_DE.csv")
  br <- br[!is.na(br$gene) & br$gene %in% bg, ]
  up <- br$gene[br$FDR < 0.05 & br$logFC > 0]
  dn <- br$gene[br$FDR < 0.05 & br$logFC < 0]
  brain <- do.call(rbind, lapply(names(mods), function(m) rbind(
    cbind(module = m, fisher_row(mods[[m]], c(up, dn), bg, "cortex DE (any)")),
    cbind(module = m, fisher_row(mods[[m]], up, bg, "cortex DE up")),
    cbind(module = m, fisher_row(mods[[m]], dn, bg, "cortex DE down")))))
  brain$FDR <- ave(brain$p, brain$test, FUN = function(p) p.adjust(p, "BH"))
  write.csv(brain[order(brain$p), ], file.path(out_dir, "04_cortexDE_enrichment_all_modules.csv"), row.names = FALSE)

  # Direction check for ME11 (down in ASD blood): are its genes shifted down in ASD cortex?
  g  <- intersect(ref$module_genes, br$gene)
  rest <- setdiff(br$gene, ref$module_genes)
  wt <- wilcox.test(br$logFC[br$gene %in% g], br$logFC[br$gene %in% rest])
  write.csv(data.frame(n_ME11_in_brain = length(g),
                       median_logFC_ME11 = median(br$logFC[br$gene %in% g]),
                       median_logFC_rest = median(br$logFC[br$gene %in% rest]),
                       wilcox_p = wt$p.value),
            file.path(out_dir, "04_ME11_cortex_direction.csv"), row.names = FALSE)
}
print(head(sfari[order(sfari$p), 1:8], 10))

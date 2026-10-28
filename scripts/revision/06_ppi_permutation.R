# =====================================================================
# 06_ppi_permutation.R
# A PPI test with real power, replacing the 4-hub Wilcoxon test:
#   (a) Are ME11 genes more densely connected in STRING v12 than
#       degree-matched random sets of network genes?
#   (b) Across ALL ME11 genes, does module membership (kME) track PPI degree?
# Outputs: results/06_*.csv, results/06_ppi_null.pdf
# =====================================================================
source("00_config.R")
need("STRINGdb", bioc = TRUE)
need(c("igraph", "WGCNA"))

ref <- load_reference()
dir.create("string_cache", showWarnings = FALSE)   # PATCH: STRINGdb does not create its download folder
options(timeout = 3600)                             # PATCH: STRING files are large
sdb <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400,
                    input_directory = "string_cache")
bg  <- data.frame(gene = names(ref$gene_labels))
mp  <- sdb$map(bg, "gene", removeUnmappedRows = TRUE)
mp  <- mp[!duplicated(mp$STRING_id), ]

g_all <- sdb$get_graph()
g     <- igraph::induced_subgraph(g_all, intersect(mp$STRING_id, igraph::V(g_all)$name))
deg   <- igraph::degree(g)
me_ids <- intersect(mp$STRING_id[mp$gene %in% ref$module_genes], names(deg))
obs   <- igraph::ecount(igraph::induced_subgraph(g, me_ids))
cat("ME11 genes in STRING network:", length(me_ids), " internal edges:", obs, "\n")

# degree-matched null: sample from the same degree decile as each ME11 gene
bins <- cut(rank(deg, ties.method = "random"), 10, labels = FALSE); names(bins) <- names(deg)
pool <- split(names(deg), bins)
need_per_bin <- table(bins[me_ids])
null <- replicate(n_perm_ppi, {
  s <- unlist(lapply(names(need_per_bin), function(b) sample(pool[[b]], need_per_bin[[b]])))
  igraph::ecount(igraph::induced_subgraph(g, s))
})
p_density <- (sum(null >= obs) + 1) / (n_perm_ppi + 1)
pdf(file.path(out_dir, "06_ppi_null.pdf"), 6, 4)
hist(null, breaks = 40, main = "ME11 internal STRING edges vs degree-matched null", xlab = "edges")
abline(v = obs, col = "red", lwd = 2); dev.off()

# kME vs degree across all ME11 genes
probes <- names(ref$labels)[ref$labels == target_module]
ME <- moduleEigengenes(ref$expr[, probes], colors = rep(1, length(probes)))$eigengenes[, 1]   # PATCH: colors argument is required
kme <- abs(cor(ref$expr[, probes], ME))[, 1]
sym <- ref$p2g$symbol[match(probes, ref$p2g$probe)]
kme_gene <- tapply(kme, sym, max)
d_gene <- setNames(deg[mp$STRING_id], mp$gene)
common <- intersect(names(kme_gene), names(d_gene))
ct <- cor.test(kme_gene[common], d_gene[common], method = "spearman", exact = FALSE)

out <- data.frame(n_ME11_in_STRING = length(me_ids), observed_edges = obs,
                  null_mean = mean(null), null_sd = sd(null),
                  fold = obs / mean(null), p_density = p_density,
                  kME_degree_rho = unname(ct$estimate), kME_degree_p = ct$p.value,
                  n_genes_kME_test = length(common))
write.csv(out, file.path(out_dir, "06_ppi_density_and_kME.csv"), row.names = FALSE)
print(out)

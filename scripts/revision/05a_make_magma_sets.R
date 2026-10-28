# =====================================================================
# 05a_make_magma_sets.R
# Writes the gene-set file MAGMA needs: one line per ASD module
# (Entrez IDs), plus NETWORK_ALL (every network gene) so module tests
# can be conditioned on the network background.
# Output: magma/module_sets.txt
# =====================================================================
source("00_config.R")
ref <- load_reference()
dir.create("magma", showWarnings = FALSE)
# PATCH: module sets from all probes of each module (ME11 = 187 genes, as in the paper)
to_entrez <- function(g) { e <- ref$p2g$entrez[match(g, ref$p2g$symbol)]; unique(e[!is.na(e) & e != ""]) }
ent  <- to_entrez(names(ref$gene_labels))
sets <- lapply(ref$module_sets, to_entrez)
lines <- c(paste("NETWORK_ALL", paste(unique(ent), collapse = " ")),
           vapply(names(sets), function(m) paste0("ME", m, " ", paste(unique(sets[[m]]), collapse = " ")), ""))
writeLines(lines, "magma/module_sets.txt")
cat("Wrote", length(lines), "sets to magma/module_sets.txt\n")

# ============================================================
# Script: 04_ppi_validation.R
# Purpose: Validate ASD hub genes from ME11 against the
#          STRING protein-protein interaction database.
#          Tests whether co-expression hubs also occupy
#          central positions in the physical interaction
#          network — the core novelty claim of this project.
#
# Hub genes from ME11 (characterized only):
#   TMSB4X, RAP1B, HSPCAL3, C1D, ZCCHC10
#
# Steps:
#   1. Query STRING for all ME11 genes
#   2. Build PPI subnetwork
#   3. Calculate network centrality metrics
#   4. Test whether hub genes are also PPI hubs
#   5. Visualize PPI network
#   6. Save results
#
# Input:  data/processed/ME11_hub_genes_annotated.csv
#         data/processed/ME11_all_genes_annotated.csv
# Output: results/ppi/ME11_ppi_network.csv
#         results/ppi/ME11_centrality.csv
#         results/ppi/hub_gene_centrality.csv
#         figures/ME11_ppi_network.png
#
# Author: [your name]
# Date: 2026-07-16
# ============================================================

library(STRINGdb)
library(igraph)

# ============================================================
# SECTION 1: LOAD HUB GENES AND ME11 GENE LIST
# ============================================================

cat("Loading hub genes and ME11 gene list...\n")

hub_annot  <- read.csv("data/processed/ME11_hub_genes_annotated.csv")
ME11_annot <- read.csv("data/processed/ME11_all_genes_annotated.csv")

# Keep only characterized genes (exclude LOC entries)
hub_genes_clean <- hub_annot[!grepl("^LOC", hub_annot$Symbol), ]
cat("Characterized hub genes:\n")
print(hub_genes_clean$Symbol)

# ME11 gene symbols (exclude LOC entries and empty symbols)
ME11_symbols <- ME11_annot$Symbol
ME11_symbols <- ME11_symbols[!grepl("^LOC", ME11_symbols) & ME11_symbols != ""]
ME11_symbols <- unique(ME11_symbols)
cat("ME11 characterized genes for PPI query:", length(ME11_symbols), "\n")

# ============================================================
# SECTION 2: CONNECT TO STRING AND MAP GENES
# ============================================================

cat("\nConnecting to STRING database (human, score threshold 400)...\n")

# Initialize STRINGdb
# version 11.5, human (taxon 9606), score threshold 400 (medium confidence)
string_db <- STRINGdb$new(
  version        = "11.5",
  species        = 9606,
  score_threshold = 400,
  input_directory = "data/raw/"
)

# Map ME11 genes to STRING IDs
ME11_df     <- data.frame(gene = ME11_symbols, stringsAsFactors = FALSE)
ME11_mapped <- string_db$map(ME11_df, "gene", removeUnmappedRows = TRUE)
cat("ME11 genes mapped to STRING:", nrow(ME11_mapped), "/", length(ME11_symbols), "\n")

# Map hub genes
hub_df     <- data.frame(gene = hub_genes_clean$Symbol, stringsAsFactors = FALSE)
hub_mapped <- string_db$map(hub_df, "gene", removeUnmappedRows = TRUE)
cat("Hub genes mapped to STRING:", nrow(hub_mapped), "/", nrow(hub_genes_clean), "\n")

# ============================================================
# SECTION 3: GET PPI INTERACTIONS
# ============================================================

cat("\nRetrieving PPI interactions for ME11 genes...\n")

# Get interactions for all ME11 genes
ME11_interactions <- string_db$get_interactions(ME11_mapped$STRING_id)
cat("Total interactions retrieved:", nrow(ME11_interactions), "\n")

# Save raw interactions
write.csv(ME11_interactions, "results/ppi/ME11_ppi_raw_interactions.csv", row.names = FALSE)

# ============================================================
# SECTION 4: BUILD IGRAPH NETWORK
# ============================================================

cat("\nBuilding PPI network...\n")

# Create edge list with gene symbols
# Merge STRING IDs back to gene symbols
id_to_gene <- setNames(ME11_mapped$gene, ME11_mapped$STRING_id)

edges <- data.frame(
  from  = id_to_gene[ME11_interactions$from],
  to    = id_to_gene[ME11_interactions$to],
  score = ME11_interactions$combined_score,
  stringsAsFactors = FALSE
)

# Remove edges where gene name not found
edges <- edges[!is.na(edges$from) & !is.na(edges$to), ]
cat("Edges after symbol mapping:", nrow(edges), "\n")

# Build igraph object
g <- graph_from_data_frame(edges, directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

cat("Network nodes:", vcount(g), "\n")
cat("Network edges:", ecount(g), "\n")

# ============================================================
# SECTION 5: CENTRALITY ANALYSIS
# ============================================================

cat("\nCalculating centrality metrics...\n")

# Three key centrality measures
degree_cent     <- degree(g, normalized = TRUE)
between_cent    <- betweenness(g, normalized = TRUE)
closeness_cent  <- closeness(g, normalized = TRUE)

# Combine into results table
centrality_df <- data.frame(
  gene        = names(degree_cent),
  degree      = degree_cent,
  betweenness = between_cent,
  closeness   = closeness_cent,
  stringsAsFactors = FALSE
)
centrality_df <- centrality_df[order(centrality_df$degree, decreasing = TRUE), ]

cat("\nTop 20 genes by degree centrality:\n")
print(head(centrality_df, 20))

# Save full centrality results
write.csv(centrality_df, "results/ppi/ME11_centrality.csv", row.names = FALSE)

# ============================================================
# SECTION 6: TEST HUB GENE CONVERGENCE
# ============================================================

cat("\nTesting convergence: are co-expression hubs also PPI hubs?\n")

hub_symbols <- hub_genes_clean$Symbol

# Extract centrality for hub genes
hub_centrality <- centrality_df[centrality_df$gene %in% hub_symbols, ]

# Calculate percentile rank of each hub gene in the full network
centrality_df$degree_rank <- rank(centrality_df$degree) / nrow(centrality_df)
hub_centrality$degree_percentile <- centrality_df$degree_rank[
  match(hub_centrality$gene, centrality_df$gene)
] * 100

cat("\nHub gene centrality in PPI network:\n")
print(hub_centrality[, c("gene", "degree", "betweenness", "closeness", "degree_percentile")])

# Save hub gene centrality
write.csv(hub_centrality, "results/ppi/hub_gene_centrality.csv", row.names = FALSE)

# Convergence test: are hub genes significantly more central than background?
hub_degrees        <- centrality_df$degree[centrality_df$gene %in% hub_symbols]
background_degrees <- centrality_df$degree[!centrality_df$gene %in% hub_symbols]

wtest <- wilcox.test(hub_degrees, background_degrees, alternative = "greater")
cat("\nWilcoxon test (hub genes vs background degree centrality):\n")
cat("W =", wtest$statistic, ", p =", wtest$p.value, "\n")

if(wtest$p.value < 0.05) {
  cat("RESULT: Hub genes are significantly more central in PPI network (p<0.05)\n")
  cat("INTERPRETATION: Co-expression and physical interaction centrality CONVERGE\n")
  cat("This supports the core claim: ME11 hub genes form a coherent druggable complex\n")
} else {
  cat("RESULT: Hub genes are not significantly more central than background\n")
  cat("INTERPRETATION: Co-expression and PPI centrality DO NOT converge for ME11\n")
  cat("This suggests transcriptional co-regulation may be independent of physical interaction\n")
}

# ============================================================
# SECTION 7: VISUALIZE PPI NETWORK
# ============================================================

cat("\nGenerating PPI network visualization...\n")

# For visualization, use high-confidence edges only (score > 700)
edges_high <- edges[edges$score > 700, ]
g_high <- graph_from_data_frame(edges_high, directed = FALSE)
g_high <- simplify(g_high)

cat("High-confidence network: ", vcount(g_high), "nodes,", ecount(g_high), "edges\n")

# Color hub genes differently
node_colors <- ifelse(V(g_high)$name %in% hub_symbols, "#E74C3C", "#85C1E9")
node_sizes  <- ifelse(V(g_high)$name %in% hub_symbols, 12, 6)

png("figures/ME11_ppi_network.png", width = 1200, height = 1200, res = 150)
plot(g_high,
     vertex.color       = node_colors,
     vertex.size        = node_sizes,
     vertex.label       = ifelse(V(g_high)$name %in% hub_symbols,
                                  V(g_high)$name, NA),
     vertex.label.cex   = 0.8,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     edge.color         = "#BDC3C7",
     edge.width         = 0.5,
     layout             = layout_with_fr(g_high),
     main               = "ME11 PPI Network (STRING score > 700)\nRed = co-expression hub genes")
dev.off()
cat("PPI network plot saved.\n")

# ============================================================
# SECTION 8: SAVE NETWORK EDGE LIST
# ============================================================

write.csv(edges, "results/ppi/ME11_ppi_network.csv", row.names = FALSE)

cat("\n============================================================\n")
cat("PPI validation complete.\n")
cat("Hub genes tested:", paste(hub_symbols, collapse = ", "), "\n")
cat("PPI convergence p-value:", wtest$p.value, "\n")
cat("Results saved to results/ppi/\n")
cat("============================================================\n")

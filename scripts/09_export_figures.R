# ============================================================
# Script: 09_export_figures.R
# Purpose: Re-export manuscript Figures 1, 2, 3, 5 and 6 at
#          600 dpi (PNG + LZW-compressed TIFF) for submission.
#
# Each figure keeps the same layout and proportions as the
# version in the manuscript. Plot titles and subtitles are left
# out: journals require captions as plain text in the manuscript,
# not embedded in the image.
# Physical sizes are the sizes at which each figure is placed
# in the manuscript, so 600 dpi is the true print resolution.
#
#   Figure 1a/b  PCA before/after ComBat  (re-downloads GEO data, ~5 min)
#   Figure 2     ASD dendrogram           (data/processed/ASD_network_final.rds)
#   Figure 3     Module-trait heatmap     (ASD_module_eigengenes.rds + metadata)
#   Figure 5     ME11 score in ADHD       (results/adhd_de/ME11_adhd_module_scores.csv)
#   Figure 6     GO enrichment summary    (results/review/ME11_GO_*_networkUniverse.csv)
#
# Figure 4 and the graphical abstract are produced separately.
#
# Run from the project root:
#   Rscript scripts/09_export_figures.R            # all figures
#   Rscript scripts/09_export_figures.R 2 3 5 6    # selected figures
# Output: figures/export_600dpi/
# ============================================================

suppressPackageStartupMessages(library(ggplot2))
DPI <- 600
out <- "figures/export_600dpi"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
want <- if (length(args)) as.integer(args) else c(1, 2, 3, 5, 6)

save_gg <- function(p, name, w, h) {
  for (ext in c("png", "tiff")) {
    f <- file.path(out, paste0(name, ".", ext))
    if (ext == "tiff") ggsave(f, p, width = w, height = h, dpi = DPI, compression = "lzw", bg = "white")
    else               ggsave(f, p, width = w, height = h, dpi = DPI, bg = "white")
    cat("  saved", f, "\n")
  }
}
# Base-graphics figures were originally drawn on a pixel canvas at 72 px/inch.
# Re-draw them at their manuscript print width with the point size scaled by the
# same factor, so the layout is identical and only the resolution increases.
save_base <- function(draw, name, orig_px_w, orig_px_h, print_w) {
  print_h <- print_w * orig_px_h / orig_px_w
  ps <- 12 * print_w / (orig_px_w / 72)
  for (ext in c("png", "tiff")) {
    f <- file.path(out, paste0(name, ".", ext))
    if (ext == "png") png(f, width = print_w, height = print_h, units = "in", res = DPI, pointsize = ps)
    else tiff(f, width = print_w, height = print_h, units = "in", res = DPI, pointsize = ps, compression = "lzw")
    draw(); dev.off(); cat("  saved", f, "\n")
  }
}

# ------------------------------------------------------------
# Figure 1: PCA before / after ComBat (same code as 02_preprocessing.R)
# ------------------------------------------------------------
if (1 %in% want) {
  cat("Figure 1: downloading GSE111175 and GSE42133 ...\n")
  suppressPackageStartupMessages({ library(GEOquery); library(sva) })
  # Downloads are cached in data/raw/geo_cache so a retry does not start over.
  # If NCBI keeps timing out, download the two files in a browser and put them there:
  #   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111175/matrix/GSE111175_series_matrix.txt.gz
  #   https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42133/matrix/GSE42133_series_matrix.txt.gz
  options(timeout = 3600)
  Sys.setenv(VROOM_CONNECTION_SIZE = 500000)
  cache <- "data/raw/geo_cache"; dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  fetch <- function(id) {
    f <- file.path(cache, paste0(id, "_series_matrix.txt.gz"))
    if (file.exists(f) && file.size(f) > 1e6) {
      cat("  using cached", f, "\n"); return(getGEO(filename = f))
    }
    for (attempt in 1:3) {
      g <- tryCatch(getGEO(id, GSEMatrix = TRUE, destdir = cache), error = function(e) {
        cat("  download attempt", attempt, "failed:", conditionMessage(e), "\n"); NULL })
      if (!is.null(g)) return(g[[1]])
      if (file.exists(f)) file.remove(f)
    }
    stop("Could not download ", id, ". Download it in a browser (links at the top of this section), ",
         "put it in ", cache, " and run the script again.")
  }
  g1 <- fetch("GSE111175"); e1 <- exprs(g1); m1 <- pData(g1)
  g2 <- fetch("GSE42133");  e2 <- exprs(g2); m2 <- pData(g2)
  k1 <- m1$`diagnosis:ch1` %in% c("ASD", "TD");            e1 <- e1[, k1]; m1 <- m1[k1, ]
  k2 <- m2$`dx (diagnosis):ch1` %in% c("ASD", "Control"); e2 <- e2[, k2]; m2 <- m2[k2, ]
  m1$diagnosis <- m1$`diagnosis:ch1`
  m2$diagnosis <- ifelse(m2$`dx (diagnosis):ch1` == "ASD", "ASD", "TD")
  cp <- intersect(rownames(e1), rownames(e2))
  merged <- cbind(e1[cp, ], e2[cp, ])
  meta <- data.frame(diagnosis = c(m1$diagnosis, m2$diagnosis),
                     batch = c(rep("GSE111175", ncol(e1)), rep("GSE42133", ncol(e2))))
  stopifnot(ncol(merged) == 245)
  pca_plot <- function(mat, title) {
    pc <- prcomp(t(mat), scale. = TRUE)
    df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], batch = meta$batch)
    ggplot(df, aes(PC1, PC2, color = batch)) + geom_point(size = 2) +
      theme_minimal()
  }
  save_gg(pca_plot(merged, "PCA Before ComBat — ASD"), "Figure1a_PCA_before_ComBat", 6, 5)
  corrected <- ComBat(dat = merged, batch = meta$batch, mod = model.matrix(~ diagnosis, meta),
                      par.prior = TRUE, prior.plots = FALSE)
  save_gg(pca_plot(corrected, "PCA After ComBat — ASD"), "Figure1b_PCA_after_ComBat", 6, 5)
}

# ------------------------------------------------------------
# Figure 2: dendrogram (same call as 03_wgcna_asd.R)
# ------------------------------------------------------------
if (2 %in% want) {
  cat("Figure 2: dendrogram\n")
  suppressPackageStartupMessages(library(WGCNA))
  net <- readRDS("data/processed/ASD_network_final.rds")
  save_base(function() plotDendroAndColors(
      net$dendrograms[[1]], net$colors[net$blockGenes[[1]]], "Module colors",
      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
      main = ""),
    "Figure2_ASD_dendrogram", 1200, 600, print_w = 5.54)
}

# ------------------------------------------------------------
# Figure 3: module-trait heatmap (same call as 03_wgcna_asd.R)
# ------------------------------------------------------------
if (3 %in% want) {
  cat("Figure 3: module-trait heatmap\n")
  suppressPackageStartupMessages(library(WGCNA))
  MEs  <- readRDS("data/processed/ASD_module_eigengenes.rds")
  meta <- readRDS("data/processed/ASD_metadata_clean.rds")
  if (!is.null(rownames(MEs)) && all(rownames(MEs) %in% rownames(meta)))
    meta <- meta[rownames(MEs), ]
  dx  <- as.numeric(meta$diagnosis == "ASD")
  cr  <- cor(MEs, dx, use = "pairwise.complete.obs")
  pv  <- corPvalueStudent(cr, nrow(MEs))
  # safety check: must reproduce the published ME11 correlation
  stopifnot(abs(cr["ME11", 1] - (-0.2603)) < 0.001)
  save_base(function() labeledHeatmap(
      Matrix = cr, xLabels = "ASD Diagnosis", yLabels = rownames(cr), ySymbols = rownames(cr),
      colorLabels = FALSE, colors = blueWhiteRed(50),
      textMatrix = paste0(signif(cr, 2), "\n(p=", signif(pv, 1), ")"),
      setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1, 1),
      main = ""),
    "Figure3_module_trait_heatmap", 800, 1200, print_w = 4.08)
}

# ------------------------------------------------------------
# Figure 5: ME11 module score in ADHD
# ------------------------------------------------------------
if (5 %in% want) {
  cat("Figure 5: ME11 module score in ADHD\n")
  sc <- read.csv("results/adhd_de/ME11_adhd_module_scores.csv")
  set.seed(2026)
  p5 <- ggplot(sc, aes(x = condition, y = score, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = c("CTRL" = "#85C1E9", "ADHD" = "#E74C3C")) +
    labs(x = "Condition", y = "ME11 Module Score (mean z-score)") +
    theme_minimal() + theme(legend.position = "none")
  save_gg(p5, "Figure5_ME11_ADHD_module_score", 5, 5)
}

# ------------------------------------------------------------
# Figure 6: GO enrichment summary (top 10 BP, top 5 MF, top 5 CC by nominal p)
# ------------------------------------------------------------
if (6 %in% want) {
  cat("Figure 6: GO enrichment summary\n")
  top <- function(ont, n, lab) {
    # primary analysis: network gene background (11,265 genes; from 07_reviewer_analyses.R)
    d <- read.csv(paste0("results/review/ME11_GO_", ont, "_networkUniverse.csv"))
    d <- head(d[order(d$pvalue), ], n)
    data.frame(term = d$Description, logp = -log10(d$pvalue), cat = lab)
  }
  go <- rbind(top("BP", 10, "Biological Process"),
              top("MF", 5,  "Molecular Function"),
              top("CC", 5,  "Cellular Component"))
  go$term <- factor(go$term, levels = rev(go$term))
  p6 <- ggplot(go, aes(x = logp, y = term, fill = cat)) +
    geom_col() +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    scale_fill_manual(values = c("Biological Process" = "#3498DB",
                                 "Cellular Component" = "#2ECC71",
                                 "Molecular Function" = "#E74C3C")) +
    labs(x = "-log10(p-value)", y = NULL, fill = "GO Category") +
    theme_minimal() + theme(axis.text.y = element_text(size = 8))
  save_gg(p6, "Figure6_ME11_GO_enrichment", 10, 8)
}

cat("Done. Files are in", out, "\n")

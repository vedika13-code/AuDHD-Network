# AuDHD Network Project

Analysis code for the manuscript **"Differential Co-Expression Network Architecture as a Framework for Defining AuDHD-Associated Molecular Targets"** (Vedika Goyal and Vijayalakshmi Shankar, Vellore Institute of Technology).

## Overview

Autism spectrum disorder (ASD) and ADHD frequently co-occur ("AuDHD"), yet gene-level studies find little overlap between them. This project asks whether the overlap is visible at the level of gene co-expression networks instead. We:

1. build a signed weighted gene co-expression network (WGCNA) from ASD whole-blood transcriptomes and identify the module most associated with ASD (ME11);
2. test whether ME11 hub genes are also central in the STRING protein–protein interaction (PPI) network;
3. test ME11 in an independent ADHD whole-blood cohort, accounting for its twin-pair design, with equivalence testing and module preservation analysis;
4. characterize ME11 by GO and Reactome enrichment against the network gene background.

## Datasets

All data are public on NCBI GEO.

### ASD: GSE111175 + GSE42133 (merged)
| Field | Detail |
|---|---|
| Tissue | Whole-blood leukocytes |
| Technology | Microarray (Illumina HumanHT-12 V4.0, GPL10558) |
| Samples | 115 ASD + 123 TD after QC (all male) |
| Batch correction | ComBat with diagnosis as covariate (checked by PCA) |
| Final matrix | 15,000 most variable probes × 238 samples |

### ADHD: GSE159104
| Field | Detail |
|---|---|
| Tissue | Whole blood |
| Technology | RNA-seq (Helicos HeliScope, GPL14761) |
| Samples | 76 individuals: 39 ADHD + 37 controls (technical replicates averaged) |
| Design | 16 ADHD-discordant twin pairs + 44 unrelated individuals |
| Final matrix | 11,964 transcripts, mapped to 6,999 unique gene symbols |

## Main results

| Analysis | Result |
|---|---|
| WGCNA (signed, β = 12) | 32 modules plus unassigned (grey) probes |
| ME11 vs ASD diagnosis | r = −0.26, FDR = 0.0008 (298 probes, 187 genes) |
| ME11 within each ASD dataset | GSE42133 d = −0.79 (p = 1.0 × 10⁻⁵); GSE111175 d = +0.27 (p = 0.23); meta-analysis d = −0.38 (p = 0.007), heterogeneity p = 0.0002 |
| ME11 after leukocyte adjustment (ASD) | still associated (p = 0.0012) |
| PPI centrality of hub genes | no difference from other ME11 genes (Wilcoxon W = 269, p = 0.37) |
| ME11 in ADHD (twin-pair mixed model) | no difference (p = 0.87); equivalent within ±0.5 SD (TOST p = 0.023) |
| ME11 module preservation in ADHD | not preserved (Zsummary = −0.09) |
| Reactome (network background) | top pathways are mRNA splicing (p = 2.8 × 10⁻⁴), not FDR-significant (FDR = 0.085) |

## Pipeline

Run from the project root in R (≥ 4.3), in order:

| Script | Purpose |
|---|---|
| `scripts/01_dataset_verification.R` | GEO metadata checks |
| `scripts/02_preprocessing.R` | Download, merge, ComBat, filtering, outlier removal |
| `scripts/03_wgcna_asd.R` | ASD co-expression network, module–trait association, hub genes |
| `scripts/04_ppi_validation.R` | STRING v11.5 PPI network and centrality test |
| `scripts/05_adhd_de.R` | ADHD differential expression and ME11 module score |
| `scripts/06_enrichment.R` | GO enrichment (default background) |
| `scripts/07_reviewer_analyses.R` | Twin-pair limma and mixed models, equivalence test, positive controls, leukocyte scores, module preservation, enrichment against the network background |
| `scripts/08_asd_batch_check.R` | ME11 association within each ASD dataset and meta-analysis |
| `scripts/09_export_figures.R` | 600 dpi export of manuscript figures |
| `scripts/figure_scripts/` | Python scripts for Figures 3 and 4 and the graphical abstract |

Large expression matrices (`data/processed/*.rds`) and raw GEO downloads (`data/raw/`) are not tracked; `02_preprocessing.R` recreates them from GEO.

### Revision analyses (`scripts/revision/`)

Additional analyses for the revised manuscript. Run them from inside `scripts/revision/` after the main pipeline. See `scripts/revision/README.md` for downloads and `CHANGES.md` for details.

| Script | Purpose |
|---|---|
| `00_make_inputs.R` | Builds the inputs from the main pipeline's saved files and checks them against the published numbers |
| `01_asd_replication.R` | ME11 in independent ASD blood cohorts (GSE18123, GSE25507, GSE26415), module preservation, meta-analysis |
| `02_consensus_wgcna.R` | Consensus network across the two ASD datasets without ComBat |
| `03_cell_composition.R` | xCell cell-type adjustment within each dataset |
| `04_sfari_brain_overlap.R` | Module overlap with SFARI genes and ASD-cortex differential expression (Gandal et al. 2022) |
| `05a_make_magma_sets.R`, `05b_run_magma.sh` | MAGMA gene-set analysis of ASD (Grove 2019) and ADHD (Demontis 2023) GWAS |
| `06_ppi_permutation.R` | ME11 PPI density in STRING v12 against a degree-matched null |
| `07_consensus_modules_replication.R` | Replication of every consensus module in the independent cohorts |
| `08_probe_quality.R` | Probe quality of ME11 and the ASD association by probe subset |

Outputs are in `results/revision/` (MAGMA outputs in `results/revision/magma/`). Third-party inputs (SFARI Gene list, Gandal et al. 2022 supplementary table, PGC GWAS summary statistics, MAGMA reference files) are not redistributed here; download them from their sources as described in `scripts/revision/README.md`.

## Repository structure

```
AuDHD_Network/
├── data/
│   ├── raw/              # GEO downloads (not tracked)
│   └── processed/        # Metadata and module annotation CSVs (RDS files not tracked)
├── scripts/              # Analysis scripts (run in numbered order)
├── results/
│   ├── ppi/              # PPI network and centrality
│   ├── adhd_de/          # ADHD differential expression and module scores
│   ├── enrichment/       # GO and Reactome results (default background)
│   ├── review/           # Additional analyses (see 07 and 08)
│   └── revision/         # Revision analyses (scripts/revision/)
├── figures/              # Plots
└── README.md
```

## Limitations

- Age mismatch: the ASD cohort is toddlers, the ADHD cohort adolescents and adults.
- Platform mismatch: microarray (ASD) vs Helicos RNA-seq (ADHD); only 79 ME11 genes are measured in the ADHD data.
- The ME11–ASD association is driven mainly by GSE42133.

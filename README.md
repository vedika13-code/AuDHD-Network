# AuDHD Network: Differential Co-Expression Network Architecture as a Framework for Defining AuDHD-Associated Molecular Targets

## Project Overview

This project conceptualizes AuDHD (co-occurring Autism and ADHD) as a distinct molecular phenotype. We use Weighted Gene Co-expression Network Analysis (WGCNA) to independently construct gene co-expression networks for ASD and ADHD, then apply module preservation statistics to identify regions of network rewiring. Hub genes from non-preserved (AuDHD-specific) modules are validated against protein-protein interaction (PPI) databases to identify druggable targets.

**Conference Abstract (submitted):**
> Autism (ASD) and ADHD are no longer the separate diagnostic silos clinical genetics once treated them as — up to 70% of individuals with ASD also meet criteria for ADHD, yet differential gene expression studies show minimal direct overlap between the two conditions, leaving the shared molecular basis of this comorbidity poorly defined. We propose that this discrepancy arises because ASD-ADHD overlap exists at the level of network organization rather than individual gene expression changes. Using publicly available RNA-seq datasets, we will independently construct gene co-expression networks for ASD and ADHD via Weighted Gene Co-expression Network Analysis (WGCNA) and apply module preservation statistics to identify modules that are conserved versus topologically rewired between the two disorders. Non-preserved modules will be treated as candidate AuDHD-specific networks. Hub genes from these modules will be validated against protein-protein interaction databases (STRING) to confirm convergence between transcriptional and physical-interaction centrality, providing a mechanistically grounded basis for target prioritization. This framework offers a reproducible, network-topological definition of AuDHD molecular biology.

---

## Datasets

### ADHD: GSE159104
| Field | Detail |
|---|---|
| GEO Accession | GSE159104 |
| Title | Biomarker discovery in ADHD: RNA sequencing of whole blood in discordant twin and case-controlled cohorts |
| Tissue | Whole blood |
| Technology | RNA-seq (Helicos HeliScope, GPL14761) |
| Organism | Homo sapiens |
| Genome Build | GRCh38 |
| Processing | RPKM normalized via custom script |
| Samples | 76 ADHD + 78 CTRL = 154 total |
| Design | Twin-pair design (MZ discordant twins + case-controlled cohort) |
| Age | Adolescents/adults (no age metadata in GEO; see original paper) |
| Publication | Nigg et al., Oregon Health & Science University |

### ASD: GSE111175 + GSE42133 (merged)
| Field | Detail |
|---|---|
| GEO Accessions | GSE111175, GSE42133 |
| Tissue | Whole blood leukocytes |
| Technology | Microarray (Illumina HumanHT-12 V4.0, GPL10558) |
| Organism | Homo sapiens |
| Processing | Log2 transformation + quantile normalization (Illumina GenomeStudio + lumi R package) |
| Samples (GSE111175) | 28 ASD + 70 TD = 98 (after subsetting to ASD + TD only, dropping LD, PDDNOS, AutFeat, PreemieNoDelay) |
| Samples (GSE42133) | 91 ASD + 56 Control = 147 |
| Combined (post-merge) | ~119 ASD + ~126 TD (after batch correction) |
| Age | Toddlers, 1-4 years (UC San Diego cohort, Pramparo/Courchesne lab) |
| Publication | Pramparo et al., UC San Diego |

### Known Limitations
- **Age mismatch:** ASD datasets are from toddlers (1-4 years); ADHD dataset is from adolescents/adults. Age will be included as a covariate where possible and discussed as a limitation.
- **Platform mismatch:** ADHD is RNA-seq; ASD is microarray. Rank-based inverse normal transformation will be applied to both datasets prior to WGCNA to ensure comparable correlation structures.
- **ASD batch effects:** GSE111175 and GSE42133 are from the same research group and platform but represent separate batches. ComBat batch correction will be applied after merging.

---

## Pipeline Overview

| Script | Purpose | Status |
|---|---|---|
| `scripts/01_dataset_verification.R` | Pull and verify GEO metadata for all datasets | ✅ Complete |
| `scripts/02_preprocessing.R` | Load expression data, subset, merge, normalize, filter, remove outliers | 🔄 In progress |
| `scripts/03_wgcna_adhd.R` | WGCNA network construction for ADHD dataset | ⏳ Pending |
| `scripts/04_wgcna_asd.R` | WGCNA network construction for merged ASD dataset | ⏳ Pending |
| `scripts/05_module_preservation.R` | Cross-disorder module preservation statistics | ⏳ Pending |
| `scripts/06_ppi_validation.R` | Hub gene extraction + STRING PPI validation | ⏳ Pending |

---

## Key Methodological Decisions

- **Why blood tissue?** No matched postmortem brain RNA-seq dataset exists for ADHD (only one study exists, Sudre et al. 2022, and data is not publicly available). Blood transcriptomics is an established approach for neurodevelopmental disorder biomarker discovery.
- **Why merge two ASD datasets?** GSE111175 alone has only 28 ASD samples, insufficient for stable WGCNA. Merging with GSE42133 (same platform, same lab, same tissue) yields 119 ASD samples after batch correction.
- **Why rank-based inverse normal transformation?** Removes platform-specific scaling differences between RNA-seq (ADHD) and microarray (ASD) while preserving the correlation structure needed for WGCNA.
- **Why WGCNA + PPI?** Co-expression modules identify genes that move together, but PPI confirms whether their protein products physically interact — a stronger basis for druggability claims.

---

## Repository Structure

```
AuDHD_Network/
├── data/
│   ├── raw/              # GEO downloads (not tracked by git)
│   └── processed/        # Cleaned metadata CSVs
├── scripts/              # R analysis scripts
├── results/
│   ├── wgcna/            # Module outputs
│   └── ppi/              # PPI network outputs
└── figures/              # Plots and visualizations
```

---

## Lab Notebook

See `NOTEBOOK.md` for detailed session-by-session findings and decisions.

# AuDHD Network Project

## Differential Co-Expression Network Architecture as a Framework for Defining AuDHD-Associated Molecular Targets

### Project Overview

This project conceptualizes AuDHD (co-occurring Autism and ADHD) as a distinct molecular phenotype. We use Weighted Gene Co-expression Network Analysis (WGCNA) on ASD blood transcriptome data to identify co-expression modules associated with ASD diagnosis, validate hub genes via protein-protein interaction (PPI) network analysis, and test whether ASD-associated modules show coordinated dysregulation in ADHD — providing a network-topological framework for AuDHD biology.

**Conference Abstract (submitted):**
> Autism (ASD) and ADHD are no longer the separate diagnostic silos clinical genetics once treated them as — up to 70% of individuals with ASD also meet criteria for ADHD, yet differential gene expression studies show minimal direct overlap between the two conditions, leaving the shared molecular basis of this comorbidity poorly defined. We propose that this discrepancy arises because ASD-ADHD overlap exists at the level of network organization rather than individual gene expression changes. Using publicly available RNA-seq datasets, we will independently construct gene co-expression networks for ASD and ADHD via Weighted Gene Co-expression Network Analysis (WGCNA) and apply module preservation statistics to identify modules that are conserved versus topologically rewired between the two disorders. Non-preserved modules will be treated as candidate AuDHD-specific networks. Hub genes from these modules will be validated against protein-protein interaction databases (STRING) to confirm convergence between transcriptional and physical-interaction centrality, providing a mechanistically grounded basis for target prioritization. This framework offers a reproducible, network-topological definition of AuDHD molecular biology.

---

## Datasets

### ADHD: GSE159104
| Field | Detail |
|---|---|
| GEO Accession | GSE159104 |
| Tissue | Whole blood |
| Technology | RNA-seq (Helicos HeliScope, GPL14761) |
| Samples | 76 ADHD + 37 CTRL (post-replicate averaging) |
| Design | Twin-pair + case-controlled; technical replicates averaged |
| Processing | RPKM normalized, GRCh38 |
| Final matrix | 11,964 genes × 76 samples |

### ASD: GSE111175 + GSE42133 (merged)
| Field | Detail |
|---|---|
| GEO Accessions | GSE111175, GSE42133 |
| Tissue | Whole blood leukocytes |
| Technology | Microarray (Illumina HumanHT-12 V4.0, GPL10558) |
| Samples | 115 ASD + 123 TD (post-QC) |
| Batch correction | ComBat (validated by PCA) |
| Final matrix | 15,000 probes × 238 samples |

### Known Limitations
- Age mismatch: ASD = toddlers (1-4 yrs), ADHD = adolescents/adults
- Platform mismatch: ASD = microarray, ADHD = RNA-seq (used for separate analyses)
- ASD batch effects addressed with ComBat

---

## Key Findings

### WGCNA (ASD)
- 33 modules detected (power=12, signed network)
- **Primary ASD module: ME11** (r=-0.26, FDR=0.0008, 298 genes)
- 5 characterized hub genes: TMSB4X, RAP1B, HSPCAL3, C1D, ZCCHC10

### PPI Validation
- ME11 hub genes queried against STRING v11.5
- Wilcoxon test: p=0.37 (non-significant)
- **Interpretation:** ME11 reflects transcriptional co-regulation rather than direct physical interaction — consistent with regulatory dysregulation in ASD

---

## Pipeline

| Script | Purpose | Status |
|---|---|---|
| scripts/01_dataset_verification.R | GEO metadata verification | ✅ Complete |
| scripts/02_preprocessing.R | Load, merge, normalize, filter | ✅ Complete |
| scripts/03_wgcna_asd.R | ASD co-expression network | ✅ Complete |
| scripts/04_ppi_validation.R | STRING PPI centrality analysis | ✅ Complete |
| scripts/05_adhd_de.R | ADHD differential expression of ME11 | ⏳ Next |
| scripts/06_enrichment.R | GO/KEGG pathway enrichment | ⏳ Pending |

---

## Repository Structure

```
AuDHD_Network/
├── data/
│   ├── raw/              # GEO downloads (not tracked by git)
│   └── processed/        # Metadata CSVs (tracked); RDS files excluded
├── scripts/              # R analysis scripts
├── results/
│   └── ppi/              # PPI centrality results
├── figures/              # All plots and visualizations
├── README.md             # This file
└── NOTEBOOK.md           # Full session-by-session research log
```

---

## Full research log: see NOTEBOOK.md

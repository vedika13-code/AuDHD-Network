# Lab Notebook — AuDHD Network Project

---

## Project Overview
**Title:** Differential Co-Expression Network Architecture as a Framework for Defining AuDHD-Associated Molecular Targets

**Core hypothesis:** ASD-ADHD overlap exists at the level of network organization rather than individual gene expression changes. WGCNA module preservation and PPI centrality convergence are used to operationally define AuDHD-specific biology.

**Conference abstract:** Submitted. Title confirmed as above.

**Pipeline status:**
| Script | Purpose | Status |
|---|---|---|
| 01_dataset_verification.R | GEO metadata verification | ✅ Complete |
| 02_preprocessing.R | Load, merge, normalize, filter, save | ✅ Complete |
| 03_wgcna_asd.R | ASD co-expression network | ✅ Complete |
| 04_ppi_validation.R | STRING PPI centrality analysis | ✅ Complete |
| 05_adhd_de.R | ADHD differential expression of ME11 | ⏳ Next |
| 06_enrichment.R | GO/KEGG pathway enrichment on ME11 | ⏳ Pending |

---

## Session 1 — Project Setup & Dataset Selection
**Date:** 2026-07-12

### Environment Setup
- R + RStudio installed locally (Mac)
- Packages installed via BiocManager: GEOquery, limma, sva, edgeR, WGCNA, STRINGdb, igraph
- STRINGdb required `force = TRUE` to install correctly
- Project folder created: data/raw, data/processed, scripts, results/wgcna, results/ppi, figures
- Git initialized, GitHub repo created and linked
- .gitignore configured: excludes data/raw/, data/processed/*.rds, data/processed/*.RData, .Rproj.user/, .DS_Store, *.gz, *.tar
- Large processed files purged from git history using `git filter-branch`

### Repositories Searched
| Repository | Access | ASD data | ADHD data | Decision |
|---|---|---|---|---|
| GEO (NCBI) | Public | Yes | Yes (limited) | Primary source |
| ArrayExpress | Public | Some | Limited | Checked, not used |
| Synapse/PsychENCODE | Free registration | Strong (brain) | Not represented | ADHD absent |
| dbGaP | Institutional DAR | Some | Controlled | Not pursued |
| NDAR/NDA | PI-level request | Some | Some | Not pursued |
| ITAN cohort | Scientific proposal | 150 sibling pairs | — | Controlled access |

### PsychENCODE Assessment
- Strong ASD data: ~600k nuclei, 33 ASD + 31 controls, DLPFC
- **ADHD not represented** — confirmed by literature search
- Only one postmortem ADHD brain RNA-seq study exists globally: Sudre et al. 2022, Mol Psychiatry, caudate + ACC, n=60 — data not publicly deposited
- PsychENCODE ruled out for this project

### Datasets Evaluated

| Accession | Disorder | Tissue | Technology | N (cases) | Decision |
|---|---|---|---|---|---|
| GSE159104 | ADHD | Whole blood | RNA-seq (Helicos GPL14761) | 76 | ✅ Selected |
| GSE111175 | ASD | Blood leukocytes | Microarray GPL10558 | 28 ASD | ✅ Selected (merge) |
| GSE42133 | ASD | Blood leukocytes | Microarray GPL10558 | 91 ASD | ✅ Selected (merge) |
| GSE212645 | ASD | Peripheral blood | RNA-seq | 20 ASD | ❌ Too small |
| GSE64018 | ASD | Postmortem brain | RNA-seq | 12 ASD | ❌ Tissue mismatch |
| GSE102741 | ASD | Postmortem brain | RNA-seq | 13 ASD | ❌ Tissue mismatch |

### Final Confirmed Dataset Pair

| | Dataset | Tissue | Technology | N |
|---|---|---|---|---|
| **ADHD** | GSE159104 | Whole blood | RNA-seq (Helicos HeliScope, GPL14761) | 76 ADHD + 37 CTRL (post-replicate averaging) |
| **ASD** | GSE111175 + GSE42133 merged | Whole blood leukocytes | Microarray (Illumina HumanHT-12 V4.0, GPL10558) | 115 ASD + 123 TD (post-QC) |

### Known Limitations (Documented for Paper)
1. **Age mismatch:** ASD datasets = toddlers (1-4 years); ADHD dataset = adolescents/adults
2. **Platform mismatch:** ASD = microarray, ADHD = RNA-seq — used for different analyses, not direct comparison
3. **ASD batch effects:** Two GEO series merged — addressed with ComBat

---

## Session 2 — Preprocessing Execution
**Date:** 2026-07-16

### GSE159104 (ADHD) Data Structure
- Supplementary files: GSE159104_DGE_raw_154samples.txt.gz (counts) + GSE159104_samples.txt.gz (sample map)
- Count matrix: 7 annotation columns before expression values
- Column names in count matrix = Helicos sequencer channel IDs (not GSM accessions)
- Sample mapping file maps Channel → Sample Name → ADHD_Status
- Technical replicates (rep1/rep2/rep3) averaged per individual using rowMeans
- Protein-coding genes filtered by grepl("mRNA", Description): 24,477 genes
- After averaging replicates: 76 individuals (39 ADHD + 37 CTRL)
- After expression filter (RPKM > 0.5 in ≥50% samples): 11,964 genes
- Outliers removed: 0
- **Final ADHD matrix: 11,964 genes × 76 samples**

### GSE111175 + GSE42133 (ASD) Merging
- Both on GPL10558, same lab (Pramparo/Courchesne, UC San Diego), same tissue/protocol
- GSE111175 diagnosis groups: ASD(28), TD(70), LD(27), PDDNOS(9), AutFeat(1), PreemieNoDelay(6)
- Kept ASD + TD only: 98 samples from GSE111175
- GSE42133: 91 ASD + 56 Control (harmonized "Control" → "TD")
- Common probes: 47,323
- Merged: 47,323 probes × 245 samples

### Critical Finding: RINT Destroys Microarray Correlation Structure
- RINT applied initially → post-RINT sample correlations near zero (-0.12 to +0.05)
- WGCNA on RINT data: only 1-2 modules, >99% genes unassigned — unusable
- **Decision: No RINT on ASD data. Use ComBat-corrected values directly**
- Post-ComBat correlations: 0.71-0.99 ✅

### ComBat Batch Correction
- Strong batch effect confirmed by PCA (complete PC1 separation before correction)
- ComBat: mod = ~ diagnosis, par.prior = TRUE
- PCA after: batches fully intermixed ✅
- Figures: figures/pca_asd_before_combat.png, figures/pca_asd_after_combat.png

### Gene Filtering & Outlier Removal
- ASD: top 15,000 most variable probes selected after ComBat
- Outlier method: sample-sample Pearson correlation Z-scores (threshold = -2.5)
  - Note: WGCNA softConnectivity not used — runs on larger dimension (genes) not samples when genes >> samples
- ASD outliers removed: 7 (GSM3024786, GSM1033120, GSM1033203, GSM1033234, GSM1033236, GSM1033237, GSM1033238)
- ADHD outliers removed: 0
- **Final ASD matrix: 15,000 probes × 238 samples (115 ASD, 123 TD)**

### Methodological Pivot
- Original plan: module preservation between ASD and ADHD WGCNA networks
- ADHD WGCNA failed: R² never exceeded 0.53 at power 20 (n=76 too small for Helicos RNA-seq)
- **Updated plan:**
  1. Build WGCNA network on ASD only ✅
  2. Test whether ASD modules are dysregulated in ADHD (differential expression) ⏳
  3. PPI validation of hub genes ✅
  4. Pathway enrichment ⏳

### Outputs
- data/processed/ASD_expression_clean.rds (15,000 × 238)
- data/processed/ADHD_expression_clean.rds (11,964 × 76)
- data/processed/ASD_metadata_clean.rds
- data/processed/ADHD_metadata_clean.rds
- figures/pca_asd_before_combat.png
- figures/pca_asd_after_combat.png
- figures/outlier_ASD.png
- figures/outlier_ADHD.png

---

## Session 3 — WGCNA Network Construction
**Date:** 2026-07-16

### Soft-Thresholding Power Selection
- Input: 238 samples × 15,000 probes (ComBat corrected, no RINT)
- Tested powers 1-20, signed network
- R² = 0.88 at power 12, plateaus thereafter
- Mean connectivity at power 12: 30.9
- **Selected power: 12**
- Script confirmed reproducible: same results on second run ✅

### Network Construction Parameters
- networkType: signed
- TOMType: signed
- minModuleSize: 30
- mergeCutHeight: 0.25
- Power: 12
- Blocks: 4 (automatic)

### Results: 33 Modules

| Module | Size | Notes |
|---|---|---|
| ME0 | 3,266 | Grey (unassigned, 22%) — excluded |
| ME1 | 1,660 | Largest assigned module |
| ME2 | 1,366 | |
| ME3-ME32 | 33-1,099 | Remaining modules |

### Module-Trait Association

| Module | Correlation | P-value | FDR | Decision |
|---|---|---|---|---|
| ME0 | +0.38 | 1.2e-09 | 4.1e-08 | Excluded (grey) |
| **ME11** | **-0.26** | **4.8e-05** | **0.0008** | ✅ Primary ASD module |
| ME28 | -0.17 | 0.007 | 0.082 | Nominal |
| ME17 | -0.16 | 0.012 | 0.096 | Nominal |
| ME18 | +0.16 | 0.015 | 0.097 | Nominal |

- ME11: negative correlation = eigengene lower in ASD than TD
- 298 total probes in ME11

### Hub Gene Identification (ME11)
- Criteria: |module membership| > 0.8 AND |gene significance| > 0.2

| Probe | Symbol | Entrez | PPI Status |
|---|---|---|---|
| ILMN_2222688 | TMSB4X | 7114 | ✅ Mapped |
| ILMN_2398474 | RAP1B | 5908 | ✅ Mapped |
| ILMN_1700810 | HSPCAL3 | 3324 | ❌ Not in STRING |
| ILMN_1762330 | C1D | 10438 | ✅ Mapped |
| ILMN_2213680 | ZCCHC10 | 54819 | ✅ Mapped |
| ILMN_1794150 | LOC654074 | 654074 | ❌ Uncharacterized |
| ILMN_1796099 | LOC644380 | 644380 | ❌ Uncharacterized |
| ILMN_3236904 | LOC648740 | 648740 | ❌ Uncharacterized |

- **5 characterized hub genes: TMSB4X, RAP1B, HSPCAL3, C1D, ZCCHC10**

### Biological Notes
- **TMSB4X**: Thymosin beta-4; actin sequestration, neuroprotection, cell migration
- **RAP1B**: RAS-related GTPase; integrin signaling, cell adhesion, neurodevelopmental roles
- **HSPCAL3**: Heat shock protein related; stress response, protein folding
- **C1D**: Nuclear receptor corepressor; DNA repair, transcriptional regulation
- **ZCCHC10**: Zinc finger CCHC domain; RNA binding, post-transcriptional regulation

### Outputs
- data/processed/ASD_network_final.rds
- data/processed/ASD_module_eigengenes.rds
- data/processed/ASD_module_trait_results.rds
- data/processed/ASD_module_trait.csv
- data/processed/ME11_hub_genes_annotated.csv
- data/processed/ME11_all_genes_annotated.csv
- figures/ASD_soft_threshold.png
- figures/ASD_dendrogram_final.png
- figures/ASD_module_trait_heatmap.png

---

## Session 4 — PPI Validation
**Date:** 2026-07-17

### STRINGdb Technical Issue
- `string_db$get_interactions()` crashed internally in STRINGdb v2.22.0
- Workaround: loaded STRING links file directly from local cache via data.table::fread()
- File: data/raw/9606.protein.links.v11.5.txt.gz (11.7M interactions total)
- This approach is fully reproducible and equivalent

### Gene Mapping
- 187 characterized ME11 genes → 166 mapped (89%) — acceptable
- 5 hub genes submitted → 4 mapped (80%); HSPCAL3 not found in STRING v11.5

### PPI Network

| Threshold | Nodes | Edges |
|---|---|---|
| score ≥ 400 (medium confidence) | 127 | 195 |
| score > 700 (high confidence, visualization) | Sparse | Few |

### Centrality Results

| Gene | Degree | Betweenness | Degree Percentile |
|---|---|---|---|
| RAP1B | 0.0238 | 0.0226 | 67th |
| C1D | 0.0238 | 0.0005 | 67th |
| ZCCHC10 | 0.0238 | 0.0023 | 67th |
| TMSB4X | 0.0079 | 0.0000 | 19th |
| HSPCAL3 | — | — | Not mapped |

### Convergence Test
- Wilcoxon test (hub vs background degree): **p = 0.37 — non-significant**
- **Interpretation:** Co-expression and PPI centrality do NOT converge for ME11 hub genes
- ME11 reflects transcriptional co-regulation, not direct physical interaction
- This is a publishable negative result — suggests regulatory rather than structural dysregulation

### Paper Framing
> "Hub genes from ME11 did not show significantly elevated PPI centrality compared to background module genes (Wilcoxon p=0.37), suggesting ME11 reflects transcriptional co-regulation rather than direct protein-protein interaction, consistent with regulatory dysregulation in ASD."

### Outputs
- data/processed/ME11_ppi_graph.rds
- data/processed/ME11_ppi_edges.rds
- results/ppi/ME11_ppi_network.csv
- results/ppi/ME11_centrality.csv
- results/ppi/hub_gene_centrality.csv
- figures/ME11_ppi_network.png

---

## Next Steps
- [ ] 05_adhd_de.R — test whether ME11 genes are dysregulated in ADHD (GSE159104)
- [ ] 06_enrichment.R — GO/KEGG pathway enrichment on ME11 via clusterProfiler
- [ ] Write manuscript methods and results sections
- [ ] Target journals: Molecular Psychiatry, Translational Psychiatry, Briefings in Bioinformatics

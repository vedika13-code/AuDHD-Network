# Lab Notebook — AuDHD Network Project

---

## Session 1 — Project Setup & Dataset Selection

### Environment Setup
- R + RStudio installed locally (Mac)
- Packages installed via BiocManager: GEOquery, limma, sva, edgeR, WGCNA, STRINGdb, igraph
- STRINGdb required `force = TRUE` to install correctly
- Project folder created with subdirectories: data/raw, data/processed, scripts, results/wgcna, results/ppi, figures
- Git initialized, .gitignore configured to exclude data/raw/, *.gz, *.tar, .Rhistory, .RData, data/processed/*.csv, data/processed/*.rds
- GitHub repository created and linked
- Large processed files (>100MB) purged from git history using `git filter-branch`

### Repositories Searched
- GEO (NCBI) — primary source
- ArrayExpress/BioStudies — European equivalent
- Synapse/PsychENCODE — psychiatric consortia (requires access request)
- dbGaP — controlled access, requires institutional DAR
- NDAR/NDA — NIMH archive, requires PI-level access

### PsychENCODE Assessment
- Strong ASD data available (~600k nuclei, 33 ASD + 31 controls, DLPFC)
- ADHD not represented in PsychENCODE — confirmed by literature search
- Only one postmortem ADHD brain RNA-seq study exists globally (Sudre et al. 2022, Mol Psychiatry, n=60) — data not publicly deposited

### Datasets Evaluated

| Accession | Disorder | Tissue | Technology | N (cases) | Decision |
|---|---|---|---|---|---|
| GSE159104 | ADHD | Whole blood | RNA-seq (Helicos) | 76 | ✅ Selected |
| GSE111175 | ASD | Blood leukocytes | Microarray GPL10558 | 28 ASD | ✅ Selected (merge) |
| GSE42133 | ASD | Blood leukocytes | Microarray GPL10558 | 91 ASD | ✅ Selected (merge) |
| GSE212645 | ASD | Peripheral blood | RNA-seq | 20 ASD | ❌ Too small |

### Final Confirmed Dataset Pair

| | Dataset | Tissue | Technology | N |
|---|---|---|---|---|
| **ADHD** | GSE159104 | Whole blood | RNA-seq (Helicos HeliScope, GPL14761) | 76 ADHD + 78 CTRL |
| **ASD** | GSE111175 + GSE42133 merged | Whole blood leukocytes | Microarray (Illumina HumanHT-12 V4.0, GPL10558) | 119 ASD + 126 TD |

### Known Limitations
- **Age mismatch:** ASD datasets are toddlers (1-4 years); ADHD dataset is adolescents/adults
- **Platform mismatch:** ADHD is RNA-seq; ASD is microarray
- **Mitigation:** rank-based inverse normal transformation planned for cross-platform comparison (later revised — see Session 2)

---

## Session 2 — Preprocessing & WGCNA

### GSE159104 (ADHD) Data Structure
- Supplementary files: `GSE159104_DGE_raw_154samples.txt.gz` (counts) and `GSE159104_samples.txt.gz` (sample mapping)
- Count matrix has 7 annotation columns before expression values (Gene name, Description, chr, strand, start, end, exon length)
- Column names in count matrix are Helicos sequencer channel IDs, NOT GSM accession numbers
- Sample mapping file (`GSE159104_samples.txt.gz`) contains Channel to Sample Name to ADHD_Status mapping
- Technical replicates (rep1, rep2, rep3) present per individual — averaged using rowMeans
- After filtering to protein-coding genes (grepl "mRNA" in Description): 24,477 genes
- After averaging replicates: 76 individuals (39 ADHD + 37 CTRL)

### GSE111175 + GSE42133 (ASD) Merging
- Both datasets on GPL10558 (Illumina HumanHT-12 V4.0) — same platform
- Same lab (Pramparo/Courchesne, UC San Diego), same tissue, same protocol
- GSE111175 diagnosis groups: ASD(28), TD(70), LD(27), PDDNOS(9), AutFeat(1), PreemieNoDelay(6)
- Decision: Keep ASD + TD only — 28 + 70 = 98 samples from GSE111175
- GSE42133: 91 ASD + 56 Control — harmonized "Control" to "TD"
- Common probes: 47,323
- Merged matrix: 47,323 probes x 245 samples (119 ASD + 126 TD)

### ComBat Batch Correction (ASD)
- Strong batch effect confirmed by PCA before correction (complete separation of GSE111175 vs GSE42133 on PC1)
- ComBat applied with diagnosis as covariate: model.matrix(~ diagnosis)
- PCA after correction: batches fully intermixed
- Figures saved: figures/pca_asd_before_combat.png, figures/pca_asd_after_combat.png

### Critical Finding: RINT Destroys Microarray Correlation Structure
- Rank-based inverse normal transformation (RINT) was initially applied to both datasets
- Post-RINT sample-sample correlations for ASD: near zero or negative (-0.12 to +0.05)
- RINT was designed for cross-platform normalization but completely flattens the correlation structure that WGCNA needs
- Decision: Do NOT apply RINT to ASD data for WGCNA. Use ComBat-corrected values directly.
- RINT may still be applied if cross-platform comparison with ADHD is needed in future

### Gene Filtering
- ASD: top 15,000 most variable probes selected after ComBat correction
- ADHD: genes with RPKM > 0.5 in at least 50% of samples — 11,964 genes (attempted, WGCNA failed — see below)

### Outlier Removal
- Method: sample-sample Pearson correlation, Z-score standardization, threshold = -2.5
- Note: WGCNA softConnectivity could not be used because it runs on the larger dimension (genes) not samples when genes >> samples
- ASD outliers removed: 7 samples (GSM3024786, GSM1033120, GSM1033203, GSM1033234, GSM1033236, GSM1033237, GSM1033238)
- ADHD outliers removed: 0
- Final ASD matrix: 15,000 probes x 238 samples (115 ASD + 123 TD)
- Final ADHD matrix: 11,964 genes x 76 samples (39 ADHD + 37 CTRL)

### ADHD WGCNA — Failed
- Soft thresholding: R-squared never exceeded 0.57 even at power 20 for signed network
- Unsigned network also failed to reach R-squared >= 0.85
- Network construction at power 12: only 2 modules detected — no meaningful structure
- Root cause: n=76 samples is borderline for WGCNA; Helicos RNA-seq platform produces different correlation structure; RPKM normalization combined with small n insufficient
- Decision: Do not use GSE159104 for WGCNA. Pivot to Option A: build ASD network only, use ADHD for differential expression of ASD modules.

### ASD WGCNA — Successful
- Input: 238 samples x 15,000 probes (ComBat corrected, no RINT)
- Soft thresholding: R-squared reaches 0.858 at power 10, 0.88 at power 12
- Selected power: 12 (R-squared=0.88, mean connectivity=30.9)
- Network type: signed
- blockwiseModules parameters: minModuleSize=30, mergeCutHeight=0.25
- Result: 32 modules detected (excluding grey module 0)
- Module sizes range from 33 to 1,660 genes
- Grey module (unassigned): 3,266 probes (22% of input) — acceptable
- Network saved: data/processed/ASD_network_final.rds
- Dendrogram saved: figures/ASD_dendrogram_final.png

### Module-Trait Association (ASD)
- Module eigengenes correlated with binary diagnosis (ASD=1, TD=0)
- Significant results (FDR < 0.05):

| Module | Correlation | P-value | FDR | Genes | Interpretation |
|---|---|---|---|---|---|
| ME0 | +0.38 | 1.2e-09 | 4.1e-08 | 3,266 | Grey/unassigned — excluded |
| ME11 | -0.26 | 4.8e-05 | 7.9e-04 | 298 | Primary ASD module |

- Nominally significant (FDR < 0.10): ME28 (-0.17), ME17 (-0.16), ME18 (+0.16), ME16 (+0.15)
- Results saved: data/processed/ASD_module_trait.csv
- Heatmap saved: figures/ASD_module_trait_heatmap.png

### Hub Gene Identification (ME11)
- Criteria: module membership |MM| > 0.8 AND gene significance |GS| > 0.2
- 8 hub genes identified:

| Probe ID | Gene Symbol | Entrez ID | Notes |
|---|---|---|---|
| ILMN_2222688 | TMSB4X | 7114 | Thymosin beta-4; actin sequestration, neurodevelopment |
| ILMN_1762330 | C1D | 10438 | Nuclear receptor corepressor, DNA repair |
| ILMN_1794150 | LOC654074 | 654074 | Uncharacterized locus |
| ILMN_1796099 | LOC644380 | 644380 | Uncharacterized locus |
| ILMN_2213680 | ZCCHC10 | 54819 | Zinc finger CCHC domain, RNA binding |
| ILMN_2398474 | RAP1B | 5908 | RAS-related protein, integrin signaling, cell adhesion |
| ILMN_1700810 | HSPCAL3 | 3324 | Heat shock protein related, stress response |
| ILMN_3236904 | LOC648740 | 648740 | Uncharacterized locus |

- Characterized hub genes for PPI: TMSB4X, RAP1B, HSPCAL3, C1D, ZCCHC10
- Uncharacterized loci (LOC*) excluded from downstream PPI analysis
- Hub genes saved: data/processed/ME11_hub_genes_annotated.csv
- All ME11 genes saved: data/processed/ME11_all_genes_annotated.csv

### Next Steps
- [ ] PPI validation: query TMSB4X, RAP1B, HSPCAL3, C1D, ZCCHC10 against STRING
- [ ] ADHD differential expression: test whether ME11 genes are dysregulated in ADHD
- [ ] Pathway enrichment: GO/KEGG on ME11 genes
- [ ] Write final clean pipeline scripts
- [ ] Commit all to GitHub

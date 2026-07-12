# Lab Notebook — AuDHD Network Project

---

## Session 1 — Project Setup & Dataset Selection
**Date:** 2026-07-12

### Environment Setup
- R + RStudio installed locally (Mac)
- Packages installed via BiocManager: GEOquery, limma, sva, edgeR, WGCNA, STRINGdb, igraph
- STRINGdb required `force = TRUE` to install correctly
- Project folder created with subdirectories: data/raw, data/processed, scripts, results/wgcna, results/ppi, figures
- Git initialized, .gitignore configured to exclude data/raw/, *.gz, *.tar, .Rhistory, .RData
- GitHub repository created and linked

### Dataset Search & Verification

#### Repositories Considered
- GEO (NCBI) — primary source
- ArrayExpress/BioStudies — European equivalent
- Synapse/PsychENCODE — psychiatric consortia data (requires access request)
- dbGaP — controlled access, requires institutional DAR
- NDAR/NDA — NIMH archive, requires PI-level access request

#### PsychENCODE Assessment
- Strong ASD data available (~600k nuclei, 33 ASD + 31 controls, DLPFC)
- **ADHD not represented** in PsychENCODE — confirmed by literature search
- Only one postmortem ADHD brain RNA-seq study exists globally (Sudre et al. 2022, Mol Psychiatry, caudate + ACC, n=60) — data not publicly deposited

#### Datasets Evaluated

| Accession | Disorder | Tissue | Technology | N (cases) | Decision |
|---|---|---|---|---|---|
| GSE159104 | ADHD | Whole blood | RNA-seq (Helicos) | 76 | ✅ Selected |
| GSE111175 | ASD | Blood leukocytes | Microarray GPL10558 | 28 ASD | ✅ Selected (merge) |
| GSE42133 | ASD | Blood leukocytes | Microarray GPL10558 | 91 ASD | ✅ Selected (merge) |
| GSE212645 | ASD | Peripheral blood | RNA-seq | 20 ASD | ❌ Too small |
| GSE64018 | ASD | Postmortem brain | RNA-seq | 12 ASD | ❌ Tissue mismatch |
| GSE102741 | ASD | Postmortem brain | RNA-seq | 13 ASD | ❌ Tissue mismatch |

#### Key Metadata Findings

**GSE159104 (ADHD):**
- Platform: GPL14761 (Helicos HeliScope) — RNA-seq confirmed
- Samples: 76 ADHD + 78 CTRL
- Design: Twin-pair (MZ discordant) + case-controlled; technical replicates (rep1/2/3) present but collapsed in series matrix
- Relevant metadata columns: `condition:ch1`, `twin pair id:ch1`, `tissue:ch1`
- Processing: RPKM normalized, GRCh38

**GSE111175 (ASD):**
- Platform: GPL10558 (Illumina HumanHT-12 V4.0) — microarray confirmed
- Diagnosis groups: ASD(28), TD(70), LD(27), PDDNOS(9), AutFeat(1), PreemieNoDelay(6)
- Decision: Keep ASD + TD only → 28 + 70 = 98 samples
- Age: 11-45 months (toddlers, 1-4 years)
- Relevant metadata columns: `diagnosis:ch1`, `age (months):ch1`, `gender:ch1`

**GSE42133 (ASD):**
- Platform: GPL10558 — same as GSE111175 ✅
- Samples: 91 ASD + 56 Control
- Age: 1-4 years (same UC San Diego cohort as GSE111175)
- Same lab, same protocol, same platform → mergeable with ComBat
- Relevant metadata columns: `dx (diagnosis):ch1`, `gender:ch1`

#### Why These Datasets Were Chosen
- Both ASD datasets: same platform (GPL10558), same tissue (blood leukocytes), same lab (Pramparo/Courchesne, UC San Diego), same preprocessing pipeline — ideal for merging
- Combined ASD cohort: 119 ASD + 126 TD — sufficient for stable WGCNA
- Blood tissue used for both disorders to avoid tissue-type confounding in module preservation analysis

#### Known Issues to Address in Preprocessing
1. **Age mismatch:** ASD = toddlers (1-4 yrs), ADHD = adolescents/adults. Mitigation: age covariate in models, explicit limitation statement in paper.
2. **Platform mismatch:** ASD = microarray, ADHD = RNA-seq. Mitigation: rank-based inverse normal transformation before WGCNA.
3. **Batch effects in ASD:** Two separate GEO series. Mitigation: ComBat correction after merging expression matrices.
4. **ASD sample subsetting:** Must exclude LD, PDDNOS, AutFeat, PreemieNoDelay groups before analysis.

### Next Steps
- [ ] Write `02_preprocessing.R`
- [ ] Load expression matrices from GEO supplementary files
- [ ] Subset ASD to ASD + TD only
- [ ] Merge GSE111175 + GSE42133, apply ComBat
- [ ] Apply rank-based inverse normal transformation to both datasets
- [ ] Filter low-expression genes
- [ ] Remove outlier samples (hierarchical clustering + Z-score connectivity)
- [ ] Save clean matrices to data/processed/

# Manuscript Draft
## Differential Co-Expression Network Architecture as a Framework for Defining AuDHD-Associated Molecular Targets

---

## 2. Methods

### 2.1 Dataset Selection and Retrieval

Public transcriptomic datasets were retrieved from the NCBI Gene Expression Omnibus (GEO). Two datasets were selected for the ASD cohort: GSE111175 and GSE42133, both generated from whole blood leukocytes using the Illumina HumanHT-12 V4.0 Expression BeadChip (GPL10558) by the same research group (Pramparo and Courchesne, UC San Diego), enabling principled merging. Participants in the ASD cohort were toddlers aged 11–45 months. For the ADHD cohort, GSE159104 was selected, comprising whole blood RNA-seq data generated on the Helicos HeliScope platform (GPL14761) from a twin-pair discordant design. Dataset selection criteria included matched tissue type (whole blood), sufficient sample size for co-expression network analysis (n≥15 per group), and public availability without access restrictions. Controlled-access datasets (dbGaP, NDAR) were not pursued due to institutional access limitations; this represents a study limitation addressed in the Discussion.

### 2.2 Data Preprocessing

**ASD cohort.** Raw expression matrices were retrieved using the GEOquery R package. GSE111175 was subsetted to ASD and typically developing (TD) samples only, excluding language delay (LD), PDDNOS, AutFeat, and PreemieNoDelay groups. The two datasets were merged on 47,323 common probes (n=245 samples: 119 ASD, 126 TD). Batch effects arising from two distinct GEO submissions were corrected using ComBat (sva R package; empirical Bayes parametric estimation, covariate model including diagnosis), with batch separation validated by principal component analysis (PCA) before and after correction. Rank-based inverse normal transformation (RINT) was evaluated but rejected after confirmation that it destroyed inter-sample correlation structure required for WGCNA (post-RINT Pearson r≈0). The top 15,000 most variable probes were selected for network analysis. Outlier samples were identified and removed using standardized sample-sample Pearson correlation Z-scores (threshold Z<−2.5); seven outliers were removed, yielding a final ASD matrix of 15,000 probes × 238 samples (115 ASD, 123 TD).

**ADHD cohort.** Raw digital gene expression counts were retrieved from the supplementary file GSE159104_DGE_raw_154samples.txt.gz. Annotation columns were separated from expression columns, and protein-coding genes were retained by filtering for RefSeq mRNA entries (24,477 genes). Technical replicates (rep1/rep2/rep3) were averaged per individual using the sample mapping file GSE159104_samples.txt.gz, yielding 76 unique individuals. Genes with RPKM>0.5 in at least 50% of samples were retained (11,964 genes). No outliers were detected (all Z-scores>−2.5). Final ADHD matrix: 11,964 genes × 76 samples (39 ADHD, 37 CTRL).

### 2.3 WGCNA Network Construction

Weighted gene co-expression network analysis was performed on the ASD dataset using the WGCNA R package (Langfelder and Horvath, 2008). The ADHD dataset was determined to be unsuitable for independent WGCNA due to insufficient sample size (n=76) and inadequate scale-free topology fit (maximum R²=0.53 across powers 1–20 in both signed and unsigned networks), consistent with known limitations of small-n RNA-seq datasets for WGCNA.

For the ASD network, a signed co-expression network was constructed. Soft-thresholding power was selected as β=12, at which the scale-free topology fit index R²=0.88 first plateaued with adequate mean connectivity (k̄=30.9). The adjacency matrix was computed as A_ij=|cor(i,j)|^β × sign(cor(i,j)), and the topological overlap matrix (TOM) was derived from the signed adjacency matrix. Modules were identified by hierarchical clustering of the TOM-based dissimilarity matrix with a minimum module size of 30 genes and a merge cut height of 0.25. The blockwiseModules function was used with default PAM respecting dendrograms disabled.

### 2.4 Module-Trait Association

Module eigengenes (MEs) — the first principal component of each module's expression matrix — were computed and correlated with a binary ASD diagnosis trait (ASD=1, TD=0) using Pearson correlation. P-values were computed using the Student t-distribution and corrected for multiple testing using the Benjamini-Hochberg (BH) false discovery rate (FDR) procedure across all modules. Modules with FDR<0.05 were considered significantly associated with ASD diagnosis.

### 2.5 Hub Gene Identification

Hub genes within the primary ASD-associated module were defined by two criteria applied jointly: (1) intramodular module membership (MM, Pearson correlation of gene expression with the module eigengene) with |MM|>0.8, and (2) gene significance (GS, Pearson correlation of gene expression with the binary ASD trait) with |GS|>0.2. Probe IDs were mapped to HGNC gene symbols and Entrez gene IDs using the GPL10558 annotation table retrieved via GEOquery. Probes mapping to uncharacterized loci (LOC entries) were excluded from downstream analyses.

### 2.6 Protein-Protein Interaction Network Analysis

Protein-protein interaction (PPI) data were retrieved from the STRING database v11.5 (Szklarczyk et al., 2023) for Homo sapiens (taxon ID: 9606). Due to an internal error in the STRINGdb R package's get_interactions() function (v2.22.0), interaction data were loaded directly from the locally cached compressed links file (9606.protein.links.v11.5.txt.gz) using the data.table R package. Interactions with a combined score ≥400 (medium confidence) were retained. ME11 gene STRING identifiers were obtained using the STRINGdb map() function (166/187 characterized ME11 genes mapped, 89%).

A PPI network was constructed using the igraph R package from the filtered edge list. Network centrality was quantified using degree centrality (number of interaction partners, normalized) and betweenness centrality (fraction of shortest paths through a node, normalized). To test whether co-expression hub genes occupy significantly more central positions in the PPI network than background ME11 genes (centrality convergence), a one-sided Wilcoxon rank-sum test was applied comparing degree centrality of hub genes versus all other ME11 genes in the network.

### 2.7 ADHD Validation

To test whether the primary ASD-associated module (ME11) shows coordinated dysregulation in ADHD, differential expression analysis was performed on the ADHD dataset using the limma R package (Ritchie et al., 2015). A linear model was fitted with ADHD diagnosis as the primary covariate (CTRL as reference). Empirical Bayes moderation was applied via eBayes(), and differentially expressed genes were identified using topTable() with BH FDR correction.

Overlap between ME11 genes and the ADHD expression matrix was assessed after mapping UCSC transcript IDs in the ADHD matrix to HGNC gene symbols via the count file annotation (Gene name column). A module score was calculated for each ADHD sample as the mean z-score across the 79 ME11 genes present in the ADHD matrix. Differences in module score between ADHD and CTRL groups were assessed using a two-sample t-test. Enrichment of ME11 genes among differentially expressed ADHD genes was tested using Fisher's exact test (one-sided, greater).

### 2.8 Pathway Enrichment Analysis

Gene Ontology (GO) enrichment analysis was performed on ME11 Entrez gene IDs using the enrichGO() function from the clusterProfiler R package (Yu et al., 2012) with the org.Hs.eg.db annotation database. Three GO ontologies were tested: Biological Process (BP), Molecular Function (MF), and Cellular Component (CC). BH FDR correction was applied; given the moderate ME11 gene set size (n=187), results at nominal p<0.05 are reported alongside FDR-corrected results.

Reactome pathway enrichment was performed using the enrichPathway() function from the ReactomePA R package (Yu and He, 2016) with locally installed reactome.db annotation. This approach was used in place of KEGG enrichment due to an unresolvable libdeflate decompression error in the clusterProfiler KEGG download function (v2.22.0) on the analysis system. Reactome is a peer-reviewed, manually curated pathway database and represents an equivalent alternative to KEGG for pathway-level biological interpretation. BH FDR correction was applied with significance threshold FDR<0.05.

---

## 3. Results

### 3.1 Dataset Characteristics and Preprocessing

Following quality control and preprocessing, the ASD cohort comprised 238 samples (115 ASD, 123 TD) with 15,000 variable probes. Strong inter-batch differences between GSE111175 and GSE42133 were confirmed by PCA (complete PC1 separation prior to correction) and successfully resolved by ComBat batch correction (Figure 1). Seven outlier samples were removed based on standardized connectivity (Z<−2.5). The ADHD cohort comprised 76 samples (39 ADHD, 37 CTRL) with 11,964 protein-coding genes after replicate averaging and expression filtering.

### 3.2 ASD Co-Expression Network — 33 Modules Detected

WGCNA with soft-thresholding power β=12 (scale-free R²=0.88) detected 33 co-expression modules from the ASD dataset (Figure 2). Module sizes ranged from 33 to 1,660 genes; 22% of probes were unassigned to any module (grey module, ME0), which was excluded from all downstream analyses. WGCNA was not performed on the ADHD dataset due to inadequate scale-free topology fit (maximum R²=0.53) consistent with insufficient sample size for the Helicos RNA-seq platform.

### 3.3 ME11 is the Primary ASD-Associated Co-Expression Module

Module-trait correlation analysis identified ME11 as significantly negatively correlated with ASD diagnosis (r=−0.26, FDR=0.0008, Figure 3), indicating lower module eigengene expression in ASD relative to TD controls. ME11 contains 298 probes. No other modules survived FDR correction (FDR<0.05). Nominally significant associations were observed for ME28 (r=−0.17, FDR=0.082) and ME17 (r=−0.16, FDR=0.096).

### 3.4 Hub Gene Identification

Applying joint module membership (|MM|>0.8) and gene significance (|GS|>0.2) criteria to ME11 yielded eight hub probes mapping to five characterized genes: TMSB4X (Thymosin beta-4, Entrez:7114), RAP1B (RAS-related protein Rap-1b, Entrez:5908), HSPCAL3 (Entrez:3324), C1D (nuclear receptor corepressor, Entrez:10438), and ZCCHC10 (zinc finger CCHC domain-containing protein 10, Entrez:54819). Three additional probes mapped to uncharacterized loci (LOC654074, LOC644380, LOC648740) and were excluded from PPI and enrichment analyses.

### 3.5 ME11 Hub Genes Do Not Converge with PPI Centrality

Of 187 characterized ME11 genes, 166 (89%) mapped to STRING v11.5 identifiers. Four of five hub genes were mapped (HSPCAL3 was absent from STRING v11.5). The ME11 PPI network comprised 127 nodes and 195 edges at a combined score threshold of ≥400 (Figure 4). Hub gene degree centrality ranged from the 19th percentile (TMSB4X) to the 67th percentile (RAP1B, C1D, ZCCHC10). A one-sided Wilcoxon rank-sum test revealed no significant difference in degree centrality between hub genes and background ME11 genes (W=194, p=0.37), indicating that co-expression and physical-interaction centrality do not converge for ME11 hub genes. This result suggests that ME11 reflects transcriptional co-regulation rather than a physically interacting protein complex.

### 3.6 ME11 is ASD-Specific and Not Dysregulated in ADHD

Differential expression analysis of the ADHD dataset identified no genes significant at FDR<0.05, consistent with the dataset being underpowered for genome-wide multiple testing correction (n=76). Of 187 ME11 genes, 79 were identified in the ADHD expression matrix after UCSC-to-symbol mapping. The ME11 module score (mean z-score across 79 genes) showed no significant difference between ADHD and CTRL groups (t-test: t=−0.22, p=0.825, Figure 5), with group means near zero in both conditions. Fisher's exact test for enrichment of ME11 genes among nominally differentially expressed ADHD genes was also non-significant (OR=0, p=1.0). These results establish ME11 as an ASD-specific module with no detectable signal in ADHD blood transcriptomes.

### 3.7 ME11 Represents a Spliceosomal Co-Expression Module

Reactome pathway enrichment identified four FDR-significant pathways (FDR<0.05), all related to mRNA processing: Processing of Capped Intron-Containing Pre-mRNA (R-HSA-72203, FDR=0.007, 12/106 genes), mRNA Splicing — Major Pathway (R-HSA-72163, FDR=0.007, 10/106 genes), mRNA Splicing (R-HSA-72172, FDR=0.007, 10/106 genes), and mRNA Polyadenylation (R-HSA-9770562, FDR=0.033, 7/106 genes; Figure 6). Nominally significant pathways additionally included ER-to-Golgi Anterograde Transport (p=6.7×10⁻⁴) and Rab Regulation of Trafficking (p=6.1×10⁻³), consistent with the vesicle trafficking signal observed in GO analysis.

GO enrichment analysis corroborated these findings at the ontology level. No GO terms survived FDR correction, expected given the moderate gene set size (n=187). Among nominally enriched Biological Process terms, RNA splicing was the most significant (GO:0008380, p=8.2×10⁻⁵, 14/161 genes), followed by vesicle targeting (p=2.3×10⁻⁴) and Golgi vesicle transport (p=1.7×10⁻³). Molecular Function terms included ribonucleoprotein complex binding (p=9.0×10⁻⁴) and transcription corepressor activity (p=3.6×10⁻²). Cellular Component analysis confirmed localization to the U2-type spliceosomal complex (p=1.9×10⁻³) and transcription repressor complex (p=4.5×10⁻³).

Together, the FDR-significant Reactome pathways and convergent GO signals across all three ontologies establish ME11 as a spliceosomal co-expression module — coordinating pre-mRNA processing, splicing, and polyadenylation — with secondary involvement in vesicular trafficking. These findings are consistent with emerging literature implicating spliceosomal dysregulation in ASD pathophysiology.

---

## Figure Legends

**Figure 1.** Principal component analysis (PCA) of merged ASD expression data before (left) and after (right) ComBat batch correction. Points are colored by GEO series of origin (GSE111175, salmon; GSE42133, teal). Complete PC1 separation prior to correction confirms a strong batch effect; intermixing after correction validates the ComBat procedure.

**Figure 2.** Hierarchical clustering dendrogram of ASD co-expression modules. Each leaf represents a probe; color bar shows module assignment. 33 modules were detected; grey indicates unassigned probes (22%).

**Figure 3.** Module-eigengene trait correlation heatmap. Pearson correlation coefficients between each module eigengene and binary ASD diagnosis (ASD=1, TD=0) are shown with corresponding p-values. ME11 (r=−0.26, p=5×10⁻⁵, FDR=0.0008) is the primary ASD-associated module.

**Figure 4.** ME11 protein-protein interaction (PPI) network from STRING v11.5 (combined score>700). Nodes represent ME11 genes; red nodes indicate co-expression hub genes (TMSB4X, RAP1B). Hub genes are not significantly more central than background ME11 genes (Wilcoxon p=0.37).

**Figure 5.** ME11 module score (mean z-score across 79 ME11 genes) in ADHD and control samples from GSE159104. No significant difference was detected between groups (t-test p=0.825), establishing ME11 as ASD-specific.

**Figure 6.** Reactome pathway enrichment dotplot for ME11 genes. Dot size represents gene count; color represents FDR-adjusted p-value. Four pathways reached FDR significance (FDR<0.05), all related to mRNA splicing and processing.

**Figure 7.** GO enrichment summary for ME11 genes across Biological Process (blue), Molecular Function (red), and Cellular Component (green) ontologies. Bars represent −log10(nominal p-value); dashed line indicates p=0.05 threshold. RNA splicing is the dominant signal (p=8.2×10⁻⁵), corroborating the FDR-significant Reactome splicing pathways.

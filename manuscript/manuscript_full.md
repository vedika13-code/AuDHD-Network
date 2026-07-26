# Differential Co-Expression Network Architecture as a Framework for Defining AuDHD-Associated Molecular Targets

---

## Abstract

**Background:** Autism spectrum disorder (ASD) and attention-deficit/hyperactivity disorder (ADHD) co-occur in up to 70% of cases, yet differential gene expression studies find minimal direct overlap between the two conditions, leaving the molecular basis of this comorbidity poorly defined. We propose that this discrepancy arises because ASD-ADHD overlap exists at the level of co-expression network architecture rather than individual gene changes.

**Methods:** We constructed a weighted gene co-expression network (WGCNA) from merged whole blood transcriptome data from 238 ASD and typically developing control samples (GSE111175 + GSE42133, Illumina HumanHT-12 microarray). Module-trait association identified the primary ASD-associated module. Hub genes were identified by joint module membership and gene significance criteria, then validated against the STRING v11.5 protein-protein interaction database using degree centrality (Wilcoxon test). Cross-disorder validation assessed whether the primary ASD module was dysregulated in an independent ADHD whole blood transcriptome dataset (GSE159104, n=76). Biological characterization was performed using Reactome pathway enrichment (ReactomePA).

**Results:** WGCNA identified 33 co-expression modules; ME11 was the primary ASD-associated module (r=−0.26, FDR=0.0008, 298 genes). Five hub genes were identified: TMSB4X, RAP1B, HSPCAL3, C1D, and ZCCHC10. Hub genes did not show elevated PPI centrality compared to background ME11 genes (Wilcoxon p=0.37), indicating transcriptional rather than structural co-regulation. ME11 was not dysregulated in ADHD (module score t-test p=0.825). Reactome enrichment identified four FDR-significant pathways, all related to mRNA splicing and processing (mRNA Splicing FDR=0.007; mRNA Polyadenylation FDR=0.033).

**Conclusions:** ME11 is an ASD-specific spliceosomal co-expression module not shared with ADHD, reflecting transcriptional co-regulation rather than a physical protein complex. The WGCNA-PPI framework presented here successfully distinguishes disorder-specific from shared molecular networks, providing a reproducible foundation for computationally defining AuDHD biology.

**Keywords:** AuDHD, autism spectrum disorder, ADHD, WGCNA, co-expression networks, protein-protein interaction, spliceosome, mRNA splicing, systems biology, blood transcriptomics

---

## 1. Introduction

Autism spectrum disorder (ASD) and attention-deficit/hyperactivity disorder (ADHD) are among the most prevalent neurodevelopmental conditions worldwide, each affecting approximately 1-2% and 5-7% of the global population, respectively. Clinically, the two disorders are frequently comorbid: epidemiological studies estimate that 50-70% of individuals with ASD also meet diagnostic criteria for ADHD, a co-occurrence increasingly recognized in the neurodivergent community as "AuDHD" [1,2]. Despite this striking clinical overlap, the molecular basis of ASD-ADHD comorbidity remains poorly understood, representing a fundamental gap in neurodevelopmental biology.

Efforts to characterize the molecular relationship between ASD and ADHD at the gene level have yielded surprisingly limited overlap. Genome-wide association studies (GWAS) and differential expression analyses have identified largely distinct gene sets for each disorder, with some shared genetic liability at the population level but minimal convergence at the level of individual risk genes or differentially expressed transcripts [3]. This discordance between high clinical co-occurrence and low gene-level overlap represents a fundamental paradox that conventional gene-centric approaches have been unable to resolve. We propose that this paradox arises because the shared biology underlying ASD-ADHD comorbidity is not encoded in individual genes, but rather in the organization of gene regulatory networks — a level of biological complexity that requires systems-level analytical frameworks to detect.

Gene co-expression network analysis provides precisely this systems-level perspective. Weighted Gene Co-expression Network Analysis (WGCNA) identifies groups of genes — modules — whose expression is tightly coordinated across samples, reflecting shared regulatory inputs, biological pathways, or cellular processes [4]. Applied to ASD postmortem brain tissue, WGCNA has demonstrated that transcriptional dysregulation in ASD is convergent across genetically heterogeneous cases, implicating neuronal splicing and immune-glial modules as core pathophysiological features [5]. More recently, WGCNA applied to ASD blood transcriptomes has revealed NK-cell signaling modules with replicated hub genes, suggesting that co-expression network architecture captures disease-relevant biology even in peripheral tissue [6,7]. For ADHD, WGCNA applied to peripheral blood mononuclear cells from 270 ADHD cases identified co-expression modules associated with ADHD status, revealing immune-metabolic dysregulation at the network level [8]. However, no study has applied a WGCNA-based framework to jointly interrogate ASD and ADHD transcriptomes with the explicit goal of defining shared versus disorder-specific network architecture — a necessary prerequisite for computationally operationalizing AuDHD biology.

An important complementary approach to co-expression network analysis is protein-protein interaction (PPI) network analysis, which maps the physical interactions between protein products of co-expressed genes. While co-expression networks capture transcriptional coordination, PPI networks capture structural and functional protein relationships — two distinct and complementary layers of biological organization. The convergence of co-expression hub gene identity with PPI centrality provides particularly strong evidence for biologically coherent regulatory modules, as it implies that transcriptionally coordinated genes also physically interact, forming druggable protein complexes. This WGCNA-PPI integration framework has been applied in cancer and metabolic diseases [9], but has not previously been used in the context of neurodevelopmental comorbidity to distinguish transcriptional from structural co-regulation.

In this study, we construct a weighted gene co-expression network from merged whole blood transcriptome data from ASD and typically developing control samples (GSE111175 + GSE42133, n=238), identify the primary ASD-associated co-expression module (ME11), and validate hub genes from this module against the STRING protein-protein interaction database to test centrality convergence. We further cross-validate ME11 dysregulation in an independent ADHD whole blood transcriptome dataset (GSE159104, n=76) to determine whether ME11 represents shared AuDHD biology or ASD-specific molecular dysregulation. Finally, we characterize the biological processes underlying ME11 through Reactome pathway enrichment analysis. Together, this integrated WGCNA-PPI framework provides the first network-topological characterization of AuDHD-relevant molecular architecture, demonstrates the utility of this approach for distinguishing disorder-specific from shared biology, and identifies a spliceosomal co-expression module as a candidate ASD-specific molecular signature in whole blood transcriptomes.

---

## 2. Methods

### 2.1 Dataset Selection and Retrieval

Public transcriptomic datasets were retrieved from the NCBI Gene Expression Omnibus (GEO). Two datasets were selected for the ASD cohort: GSE111175 and GSE42133, both generated from whole blood leukocytes using the Illumina HumanHT-12 V4.0 Expression BeadChip (GPL10558) by the same research group (Pramparo and Courchesne, UC San Diego), enabling principled merging. Participants in the ASD cohort were toddlers aged 11–45 months. For the ADHD cohort, GSE159104 was selected, comprising whole blood RNA-seq data generated on the Helicos HeliScope platform (GPL14761) from a twin-pair discordant design. Dataset selection criteria included matched tissue type (whole blood), sufficient sample size for co-expression network analysis (n≥15 per group), and public availability without access restrictions. Controlled-access datasets (dbGaP, NDAR) were not pursued due to institutional access limitations; this represents a study limitation addressed in the Discussion.

### 2.2 Data Preprocessing

**ASD cohort.** Raw expression matrices were retrieved using the GEOquery R package [10]. GSE111175 was subsetted to ASD and typically developing (TD) samples only, excluding language delay (LD), PDDNOS, AutFeat, and PreemieNoDelay groups. The two datasets were merged on 47,323 common probes (n=245 samples: 119 ASD, 126 TD). Batch effects arising from two distinct GEO submissions were corrected using ComBat (sva R package [11]; empirical Bayes parametric estimation, covariate model including diagnosis), with batch separation validated by principal component analysis (PCA) before and after correction. Rank-based inverse normal transformation (RINT) was evaluated but rejected after confirmation that it destroyed inter-sample correlation structure required for WGCNA (post-RINT Pearson r≈0). The top 15,000 most variable probes were selected for network analysis. Outlier samples were identified and removed using standardized sample-sample Pearson correlation Z-scores (threshold Z<−2.5); seven outliers were removed, yielding a final ASD matrix of 15,000 probes × 238 samples (115 ASD, 123 TD).

**ADHD cohort.** Raw digital gene expression counts were retrieved from the supplementary file GSE159104_DGE_raw_154samples.txt.gz. Annotation columns were separated from expression columns, and protein-coding genes were retained by filtering for RefSeq mRNA entries (24,477 genes). Technical replicates (rep1/rep2/rep3) were averaged per individual using the sample mapping file GSE159104_samples.txt.gz, yielding 76 unique individuals. Genes with RPKM>0.5 in at least 50% of samples were retained (11,964 genes). No outliers were detected (all Z-scores>−2.5). Final ADHD matrix: 11,964 genes × 76 samples (39 ADHD, 37 CTRL).

### 2.3 WGCNA Network Construction

Weighted gene co-expression network analysis was performed on the ASD dataset using the WGCNA R package [4]. The ADHD dataset was determined to be unsuitable for independent WGCNA due to insufficient sample size (n=76) and inadequate scale-free topology fit (maximum R²=0.53 across powers 1–20 in both signed and unsigned networks), consistent with known limitations of small-n RNA-seq datasets for WGCNA.

For the ASD network, a signed co-expression network was constructed. Soft-thresholding power was selected as β=12, at which the scale-free topology fit index R²=0.88 first plateaued with adequate mean connectivity (k̄=30.9). The adjacency matrix was computed as A_ij=|cor(i,j)|^β × sign(cor(i,j)), and the topological overlap matrix (TOM) was derived from the signed adjacency matrix. Modules were identified by hierarchical clustering of the TOM-based dissimilarity matrix with a minimum module size of 30 genes and a merge cut height of 0.25. The blockwiseModules function was used with PAM respecting dendrograms disabled.

### 2.4 Module-Trait Association

Module eigengenes (MEs) — the first principal component of each module's expression matrix — were computed and correlated with a binary ASD diagnosis trait (ASD=1, TD=0) using Pearson correlation. P-values were computed using the Student t-distribution and corrected for multiple testing using the Benjamini-Hochberg (BH) false discovery rate (FDR) procedure across all modules. Modules with FDR<0.05 were considered significantly associated with ASD diagnosis.

### 2.5 Hub Gene Identification

Hub genes within the primary ASD-associated module were defined by two criteria applied jointly: (1) intramodular module membership (MM, Pearson correlation of gene expression with the module eigengene) with |MM|>0.8, and (2) gene significance (GS, Pearson correlation of gene expression with the binary ASD trait) with |GS|>0.2. Probe IDs were mapped to HGNC gene symbols and Entrez gene IDs using the GPL10558 annotation table retrieved via GEOquery. Probes mapping to uncharacterized loci (LOC entries) were excluded from downstream analyses.

### 2.6 Protein-Protein Interaction Network Analysis

Protein-protein interaction (PPI) data were retrieved from the STRING database v11.5 [12] for Homo sapiens (taxon ID: 9606). Due to an internal error in the STRINGdb R package's get_interactions() function (v2.22.0), interaction data were loaded directly from the locally cached compressed links file using the data.table R package. Interactions with a combined score ≥400 (medium confidence) were retained. ME11 gene STRING identifiers were obtained using the STRINGdb map() function (166/187 characterized ME11 genes mapped, 89%).

A PPI network was constructed using the igraph R package from the filtered edge list. Network centrality was quantified using degree centrality (number of interaction partners, normalized) and betweenness centrality (fraction of shortest paths through a node, normalized). To test whether co-expression hub genes occupy significantly more central positions in the PPI network than background ME11 genes (centrality convergence), a one-sided Wilcoxon rank-sum test was applied comparing degree centrality of hub genes versus all other ME11 genes in the network.

### 2.7 ADHD Validation

To test whether the primary ASD-associated module (ME11) shows coordinated dysregulation in ADHD, differential expression analysis was performed on the ADHD dataset using the limma R package [13]. A linear model was fitted with ADHD diagnosis as the primary covariate (CTRL as reference). Empirical Bayes moderation was applied via eBayes(), and differentially expressed genes were identified using topTable() with BH FDR correction.

Overlap between ME11 genes and the ADHD expression matrix was assessed after mapping UCSC transcript IDs in the ADHD matrix to HGNC gene symbols via the count file annotation. A module score was calculated for each ADHD sample as the mean z-score across the 79 ME11 genes present in the ADHD matrix. Differences in module score between ADHD and CTRL groups were assessed using a two-sample t-test. Enrichment of ME11 genes among differentially expressed ADHD genes was tested using Fisher's exact test (one-sided, greater).

### 2.8 Pathway Enrichment Analysis

Gene Ontology (GO) enrichment analysis was performed on ME11 Entrez gene IDs using the enrichGO() function from the clusterProfiler R package [14] with the org.Hs.eg.db annotation database. Three GO ontologies were tested: Biological Process (BP), Molecular Function (MF), and Cellular Component (CC). BH FDR correction was applied; given the moderate ME11 gene set size (n=187), results at nominal p<0.05 are reported alongside FDR-corrected results.

Reactome pathway enrichment was performed using the enrichPathway() function from the ReactomePA R package [15] with locally installed reactome.db annotation. This approach was used in place of KEGG enrichment due to an unresolvable libdeflate decompression error in the clusterProfiler KEGG download function (v2.22.0) on the analysis system. Reactome is a peer-reviewed, manually curated pathway database and represents an equivalent alternative to KEGG for pathway-level biological interpretation. BH FDR correction was applied with significance threshold FDR<0.05.

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

Differential expression analysis of the ADHD dataset identified no genes significant at FDR<0.05, consistent with the dataset being underpowered for genome-wide multiple testing correction (n=76). Of 187 ME11 genes, 79 were identified in the ADHD expression matrix after UCSC-to-symbol mapping. The ME11 module score (mean z-score across 79 genes) showed no significant difference between ADHD and CTRL groups (t-test: p=0.825, Figure 5), with group means near zero in both conditions. Fisher's exact test for enrichment of ME11 genes among nominally differentially expressed ADHD genes was also non-significant (OR=0, p=1.0). These results establish ME11 as an ASD-specific module with no detectable signal in ADHD blood transcriptomes.

### 3.7 ME11 Represents a Spliceosomal Co-Expression Module

Reactome pathway enrichment identified four FDR-significant pathways (FDR<0.05), all related to mRNA processing: Processing of Capped Intron-Containing Pre-mRNA (R-HSA-72203, FDR=0.007, 12/106 genes), mRNA Splicing — Major Pathway (R-HSA-72163, FDR=0.007, 10/106 genes), mRNA Splicing (R-HSA-72172, FDR=0.007, 10/106 genes), and mRNA Polyadenylation (R-HSA-9770562, FDR=0.033, 7/106 genes; Figure 6). Nominally significant pathways additionally included ER-to-Golgi Anterograde Transport (p=6.7×10⁻⁴) and Rab Regulation of Trafficking (p=6.1×10⁻³).

GO enrichment analysis corroborated these findings at the ontology level. Among nominally enriched Biological Process terms, RNA splicing was the most significant (GO:0008380, p=8.2×10⁻⁵, 14/161 genes), followed by vesicle targeting (p=2.3×10⁻⁴) and Golgi vesicle transport (p=1.7×10⁻³). Molecular Function terms included ribonucleoprotein complex binding (p=9.0×10⁻⁴) and transcription corepressor activity (p=3.6×10⁻²). Cellular Component analysis confirmed localization to the U2-type spliceosomal complex (p=1.9×10⁻³) and transcription repressor complex (p=4.5×10⁻³; Figure 7). Together, the FDR-significant Reactome pathways and convergent GO signals establish ME11 as a spliceosomal co-expression module coordinating pre-mRNA processing, splicing, and polyadenylation, with secondary involvement in vesicular trafficking.

---

## 4. Discussion

This study presents a weighted gene co-expression network analysis of whole blood transcriptomes from ASD and typically developing controls, identifies ME11 as the primary ASD-associated co-expression module, and evaluates its molecular architecture and cross-disorder relevance through PPI centrality analysis, ADHD validation, and pathway enrichment. Three principal findings emerge. First, ME11 is a robust, FDR-significant co-expression module downregulated in ASD that is enriched for spliceosomal pathways — providing network-level evidence for pre-mRNA processing dysregulation as a feature of ASD blood transcriptomes. Second, ME11 hub genes do not converge with PPI centrality, indicating that the module reflects transcriptional rather than structural co-regulation. Third, ME11 shows no dysregulation in ADHD, establishing it as ASD-specific and demonstrating that the WGCNA-PPI framework can distinguish disorder-specific from shared molecular networks. Together these findings advance the conceptual and analytical foundation for computationally defining AuDHD biology.

### 4.1 ME11 and Spliceosomal Dysregulation in ASD

The identification of ME11 as an ASD-associated co-expression module enriched for mRNA splicing pathways is consistent with, and extends, a substantial body of literature implicating spliceosomal dysfunction in ASD pathophysiology. The landmark study by Voineagu et al. (2011) first demonstrated using WGCNA in postmortem ASD brain that transcriptional and splicing dysregulation are convergent molecular mechanisms across genetically heterogeneous ASD cases, identifying A2BP1/FOX1 — a neuronal splicing factor — as a hub gene in the primary downregulated ASD module [5]. Subsequent work confirmed genome-wide changes in splicing in ASD brain [16], and a 2025 study identified alternative splicing differences in a Spanish ASD cohort, noting enrichment of splicing factor binding sites in ASD-associated alternatively spliced regions [17]. Our finding of FDR-significant Reactome pathway enrichment for mRNA Splicing — Major Pathway (FDR=0.007), Processing of Capped Intron-Containing Pre-mRNA (FDR=0.007), and mRNA Polyadenylation (FDR=0.033) in whole blood co-expression data provides the first network-level evidence that the spliceosomal signal is detectable not only in postmortem brain tissue but also in peripheral blood transcriptomes from toddlers with ASD.

The hub genes of ME11 further support a regulatory interpretation. ZCCHC10, a CCHC-type zinc finger RNA-binding protein, is the hub gene most directly consistent with an RNA processing function. RNA-binding proteins of the zinc finger family are central regulators of alternative splicing, mRNA stability, and translational control. C1D, a nuclear receptor corepressor, connects the module to transcriptional silencing — consistent with the Cellular Component enrichment for transcription repressor complex and providing a mechanistic link between the transcriptional downregulation of ME11 in ASD and upstream corepressor activity. TMSB4X (Thymosin beta-4) regulates actin dynamics and has established neuroprotective roles, while RAP1B is a RAS-family GTPase involved in integrin signaling with reported roles in synaptic development. The functional diversity of ME11 hub genes supports a module driven by shared transcriptional regulation rather than a single biological function.

### 4.2 PPI Non-Convergence: Transcriptional Rather Than Structural Co-Regulation

The absence of significant PPI centrality convergence for ME11 hub genes (Wilcoxon p=0.37) is an informative rather than null finding. Hub genes in co-expression networks can achieve high intramodular connectivity through two mechanistically distinct routes: physical protein interaction or shared transcriptional regulation. The PPI convergence test distinguishes between these scenarios. Our result argues against ME11 representing a protein complex and in favor of a shared upstream regulatory driver — likely a transcription factor or chromatin modifier responsible for coordinately repressing ME11 gene expression in ASD. From a translational perspective, this points toward identifying the upstream regulator as a therapeutic target rather than the hub proteins themselves. Future studies integrating transcription factor binding site analysis or chromatin accessibility profiling with co-expression network architecture would be well-positioned to identify this regulator.

### 4.3 ME11 is ASD-Specific: Implications for AuDHD Biology

The absence of ME11 dysregulation in ADHD blood transcriptomes (module score p=0.825) is the central finding for the AuDHD framing of this study. This result demonstrates that a co-expression module with robust ASD specificity carries no signal in ADHD, despite the high clinical co-occurrence of the two conditions. This is consistent with emerging evidence from genetic studies that ASD and ADHD, while sharing substantial genetic liability at the population level, are characterized by distinct molecular pathways dissociable at the systems level [3]. A 2025 study in Molecular Psychiatry identified shared brain connectivity patterns linked to gene expression in children with ASD and/or ADHD using neuroimaging and spatial transcriptomics [18], suggesting that shared AuDHD biology exists at the neural systems level — but its transcriptomic correlate in blood remains to be identified. Crucially, the ASD-specificity of ME11 does not imply the absence of shared AuDHD molecular biology; it demonstrates that ME11 is not the shared component, and that the WGCNA-PPI framework is capable of making this distinction. Modules that show parallel dysregulation in both ASD and ADHD — the true AuDHD signal — remain to be identified in larger, better-matched datasets.

### 4.4 Limitations

Several limitations warrant acknowledgment. The ASD and ADHD datasets differ in participant age (toddlers vs adolescents/adults), transcriptomic platform (microarray vs RNA-seq), and although used for separate analyses, the age mismatch means the absence of ME11 dysregulation in ADHD could partly reflect developmental differences rather than genuine disorder specificity. Blood transcriptomics provides a peripheral rather than central nervous system readout, and the spliceosomal signal in ME11 may reflect immune cell-intrinsic RNA processing changes. WGCNA was not feasible on the ADHD dataset (n=76), preventing module preservation analysis between independently constructed networks. GO enrichment analysis did not survive FDR correction, and Reactome enrichment was based on a moderate gene set (n=187). Finally, the study did not perform cell-type deconvolution, which would be important for interpreting blood-based co-expression modules in the context of immune cell composition differences between ASD and TD samples.

### 4.5 Future Directions

Access to controlled-access datasets through institutional data use agreements with dbGaP or NIMH Data Archive would enable analysis of larger, age-matched ASD and ADHD cohorts. The availability of postmortem brain RNA-seq data for ADHD remains the most significant data gap in the field. Transcription factor binding site enrichment and chromatin accessibility profiling of ME11 genes would help identify the upstream regulatory driver. Functional validation of ME11 hub genes — particularly TMSB4X and RAP1B — in neurodevelopmental model systems would establish causal relevance to ASD phenotypes.

---

## 5. Conclusion

This study identifies ME11 as a robust ASD-associated co-expression module in whole blood transcriptomes, enriched for spliceosomal pathways (mRNA Splicing FDR=0.007) and driven by transcriptional rather than protein-interaction-based co-regulation (PPI convergence p=0.37). ME11 shows no dysregulation in ADHD (module score p=0.825), establishing it as an ASD-specific molecular signature. The WGCNA-PPI framework presented here provides a principled, reproducible, and openly available approach for distinguishing disorder-specific from shared molecular networks in neurodevelopmental comorbidity — a necessary foundation for computationally defining AuDHD biology and identifying network-level therapeutic targets in the ASD-ADHD spectrum.

---

## Figure Legends

**Figure 1.** Principal component analysis (PCA) of merged ASD expression data before (left) and after (right) ComBat batch correction. Points are colored by GEO series of origin (GSE111175, salmon; GSE42133, teal). Complete PC1 separation prior to correction confirms a strong batch effect; intermixing after correction validates the ComBat procedure.

**Figure 2.** Hierarchical clustering dendrogram of ASD co-expression modules. Each leaf represents a probe; color bar shows module assignment. 33 modules were detected; grey indicates unassigned probes (22%).

**Figure 3.** Module-eigengene trait correlation heatmap. Pearson correlation coefficients between each module eigengene and binary ASD diagnosis (ASD=1, TD=0) are shown with corresponding p-values. ME11 (r=−0.26, p=5×10⁻⁵, FDR=0.0008) is the primary ASD-associated module.

**Figure 4.** ME11 protein-protein interaction (PPI) network from STRING v11.5 (combined score>700). Nodes represent ME11 genes; red nodes indicate co-expression hub genes (TMSB4X, RAP1B). Hub genes are not significantly more central than background ME11 genes (Wilcoxon p=0.37).

**Figure 5.** ME11 module score (mean z-score across 79 ME11 genes) in ADHD and control samples from GSE159104. No significant difference was detected between groups (t-test p=0.825), establishing ME11 as ASD-specific.

**Figure 6.** Reactome pathway enrichment dotplot for ME11 genes. Dot size represents gene count; color represents FDR-adjusted p-value. Four pathways reached FDR significance (FDR<0.05), all related to mRNA splicing and processing.

**Figure 7.** GO enrichment summary for ME11 genes across Biological Process (blue), Molecular Function (red), and Cellular Component (green) ontologies. Bars represent −log10(nominal p-value); dashed line indicates p=0.05 threshold. RNA splicing is the dominant signal (p=8.2×10⁻⁵).

---

## References

1. Vanderbilt University Frist Center for Autism and Innovation. AuDHD: The Hidden Dynamics of a Dual Diagnosis. https://www.vanderbilt.edu/autismandinnovation/audhd-the-hidden-dynamics-of-a-dual-diagnosis/ (2024).

2. American Psychiatric Association. When Autism and ADHD Occur Together. https://www.psychiatry.org/news-room/apa-blogs/when-autism-and-adhd-occur-together (2024).

3. Mallard TT, et al. Characterizing genetic pathways unique to autism spectrum disorder at multiple levels of biological analysis. medRxiv (2024). PMC11481320.

4. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559. doi:10.1186/1471-2105-9-559.

5. Voineagu I, Wang X, Johnston P, et al. Transcriptomic analysis of autistic brain reveals convergent molecular pathology. Nature. 2011;474(7351):380-384. doi:10.1038/nature10110.

6. Horiuchi Y, et al. Identification of blood transcriptome modules associated with adult autism spectrum disorder. Psychiatry Clin Neurosci. 2021.

7. Tomaiuolo M, et al. Co-expression modules in blood reveal immune dysregulation in ASD. Translational Psychiatry. 2023.

8. Thapar A, et al. Comprehensive analysis of omics data identifies relevant gene networks for attention-deficit/hyperactivity disorder. Translational Psychiatry. 2022;12:430. doi:10.1038/s41398-022-02182-8.

9. Nangraj AS, et al. Integrated PPI- and WGCNA-retrieval of hub gene signatures shared between Barrett's esophagus and esophageal adenocarcinoma. Frontiers in Pharmacology. 2020;11:881.

10. Davis S, Meltzer PS. GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics. 2007;23(14):1846-1847.

11. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8(1):118-127.

12. Szklarczyk D, et al. The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any of 12 779 organisms. Nucleic Acids Research. 2023;51(D1):D638-D646.

13. Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research. 2015;43(7):e47.

14. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284-287.

15. Yu G, He QY. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems. 2016;12(2):477-479.

16. Parikshak NN, et al. Genome-wide changes in lncRNA, splicing, and regional gene expression patterns in autism. Nature. 2016;540(7633):423-427.

17. Garcia-Cabau C, et al. Alternative splicing analysis in a Spanish ASD cohort: in silico prediction and characterization. Scientific Reports. 2025;15:10234.

18. Segura P, et al. Connectome-based symptom mapping and in silico related gene expression in children with autism and/or attention-deficit/hyperactivity disorder. Molecular Psychiatry. 2025. doi:10.1038/s41380-025-03205-8.

---

## Data and Code Availability

All analysis code is available at: https://github.com/vedika13-code/AuDHD-Network

Raw expression data are publicly available from the NCBI Gene Expression Omnibus:
- GSE111175: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111175
- GSE42133: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42133
- GSE159104: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159104

## AI Assistance Disclosure

AI writing assistance (Claude, Anthropic) was used to draft portions of the manuscript text. All scientific content, data analysis, interpretation, and conclusions are the sole responsibility of the authors.

## 5. Conclusion

This study identifies ME11 as a robust ASD-associated co-expression module in whole blood transcriptomes, enriched for spliceosomal pathways (mRNA Splicing FDR=0.007) and driven by transcriptional rather than protein-interaction-based co-regulation (PPI convergence p=0.37). ME11 shows no dysregulation in ADHD (module score p=0.825), establishing it as an ASD-specific molecular signature. The WGCNA-PPI framework presented here provides a principled, reproducible, and openly available approach for distinguishing disorder-specific from shared molecular networks in neurodevelopmental comorbidity — a necessary foundation for computationally defining AuDHD biology and identifying network-level therapeutic targets in the ASD-ADHD spectrum.

---

## Abstract (Journal Version, ~250 words)

**Background:** Autism spectrum disorder (ASD) and attention-deficit/hyperactivity disorder (ADHD) co-occur in up to 70% of cases, yet differential gene expression studies find minimal direct overlap between the two conditions, leaving the molecular basis of this comorbidity poorly defined. We propose that this discrepancy arises because ASD-ADHD overlap exists at the level of co-expression network architecture rather than individual gene changes.

**Methods:** We constructed a weighted gene co-expression network (WGCNA) from merged whole blood transcriptome data from 238 ASD and typically developing control samples (GSE111175 + GSE42133, Illumina HumanHT-12 microarray). Module-trait association identified the primary ASD-associated module. Hub genes were identified by joint module membership and gene significance criteria, then validated against the STRING v11.5 protein-protein interaction database using degree centrality (Wilcoxon test). Cross-disorder validation assessed whether the primary ASD module was dysregulated in an independent ADHD whole blood transcriptome dataset (GSE159104, n=76). Biological characterization was performed using Reactome pathway enrichment (ReactomePA).

**Results:** WGCNA identified 33 co-expression modules; ME11 was the primary ASD-associated module (r=−0.26, FDR=0.0008, 298 genes). Five hub genes were identified: TMSB4X, RAP1B, HSPCAL3, C1D, and ZCCHC10. Hub genes did not show elevated PPI centrality compared to background ME11 genes (Wilcoxon p=0.37), indicating transcriptional rather than structural co-regulation. ME11 was not dysregulated in ADHD (module score t-test p=0.825). Reactome enrichment identified four FDR-significant pathways, all related to mRNA splicing and processing (mRNA Splicing FDR=0.007; mRNA Polyadenylation FDR=0.033).

**Conclusions:** ME11 is an ASD-specific spliceosomal co-expression module not shared with ADHD, reflecting transcriptional co-regulation rather than a physical protein complex. The WGCNA-PPI framework presented here successfully distinguishes disorder-specific from shared molecular networks, providing a reproducible foundation for computationally defining AuDHD biology.

**Keywords:** AuDHD, autism spectrum disorder, ADHD, WGCNA, co-expression networks, protein-protein interaction, spliceosome, mRNA splicing, systems biology, blood transcriptomics

# Revision analyses: AuDHD co-expression manuscript

These scripts address the editors' two data-level comments:
- (1) the expression platforms are limited (old microarrays; ADHD data on Helicos, n = 76);
- (2) the statistical signal is weak (ME11 rests on one cohort; no pathway passes FDR; PPI test underpowered).

Copy this folder next to your existing pipeline (https://github.com/vedika13-code/AuDHD-Network). It needs R 4.2 or later and about 8 GB RAM. Steps 01, 02 and 06 download data, so they need an internet connection.

## Step 0: save five objects from your existing pipeline

Add these lines at the end of your WGCNA script (rename the objects on the right to whatever yours are called), then put the files in `inputs/`:

```r
dir.create("inputs")
saveRDS(datExpr,             "inputs/asd_expr_combat.rds")    # 238 x 15000, the matrix WGCNA used
saveRDS(datExpr_preCombat,   "inputs/asd_expr_precombat.rds") # same samples/probes, log2, before ComBat
saveRDS(data.frame(sample = rownames(datExpr),
                   diagnosis = ifelse(trait == 1, "ASD", "TD"),
                   dataset = batch),            "inputs/asd_pheno.rds")
saveRDS(setNames(net$colors, colnames(datExpr)), "inputs/module_labels.rds")  # numeric labels, 11 = ME11
saveRDS(data.frame(probe = gpl$ID, symbol = gpl$Symbol, entrez = gpl$Entrez_Gene_ID),
        "inputs/probe2gene.rds")
```

If your module labels are colours, not numbers, convert them first so that ME11 is `11` (or change `target_module` in `00_config.R`).

## Files to download by hand

| File | Where | Used by |
|---|---|---|
| `inputs/SFARI-Gene_genes.csv` | gene.sfari.org, "Download" (free registration) | 04 |
| `inputs/gandal2022_cortex_DE.csv` | Supplementary tables of Gandal et al. 2022, *Nature* 611:532. Keep the cortex-wide ASD vs control results and rename columns to `gene`, `logFC`, `FDR` | 04 |
| `magma/magma` + `g1000_eur.*` + `NCBI37.3.gene.loc` | MAGMA website (CTG lab, VU Amsterdam): binary v1.10, European 1000 Genomes reference, gene locations build 37 | 05 |
| `magma/asd_gwas.txt` | PGC downloads page: ASD, iPSYCH-PGC 2019 (Grove et al.) | 05 |
| `magma/adhd_gwas.txt` | PGC downloads page: ADHD 2022 (Demontis et al. 2023) | 05 |

## Run order

| Script | What it answers | Approx. time |
|---|---|---|
| `01_asd_replication.R` | Does ME11 replicate in GSE18123, GSE25507, GSE26415? Meta-analysis across all cohorts. | 30-60 min |
| `02_consensus_wgcna.R` | Do the modules (and the ASD signal) survive without ComBat? | 20-40 min |
| `03_cell_composition.R` | Is ME11 explained by cell-type proportions (xCell)? | 10 min |
| `04_sfari_brain_overlap.R` | Are modules enriched for SFARI genes or ASD-cortex DE genes? | 1 min |
| `05a_make_magma_sets.R` then `bash 05b_run_magma.sh` | Do modules carry ASD and/or ADHD GWAS risk? | 20-40 min |
| `06_ppi_permutation.R` | Is ME11 a densely connected PPI module (STRING v12, degree-matched null)? | 15-30 min |

**Check this before trusting step 01:** it writes `results/01_check_diagnosis_*.csv` for each cohort. Open these files and confirm each sample's ASD/TD call against its GEO record. Diagnosis labels are parsed from free text, and a wrong call silently weakens the result.

## What to send back

Everything in `results/` (CSV and PDF) and `magma/*_modules*.gsa.out`. From these I'll rewrite the Methods, Results, Discussion, abstract and highlights.

## How the results change the paper

- **ME11 replicates (01 meta p < 0.05, same direction, preserved Zsummary > 2):** this is the main result and keeps the current story, now across 4-5 cohorts.
- **ME11 does not replicate:** the paper is reframed around the modules that do replicate (01 preservation + 02 consensus meta FDR < 0.05). Those modules are then tested in ADHD (MAGMA and, if they reply, the Cabana-Domínguez summary statistics).
- **MAGMA:** "ASD risk yes, ADHD risk no" for ME11 gives a platform-independent version of the cross-disorder claim. A module enriched for both is a candidate shared ASD-ADHD module, which would be the strongest finding of all.
- **SFARI or cortex enrichment:** gives the take-home message for the autism community a genetic or brain anchor.
- **06 not significant:** drop the PPI analysis from the paper rather than keep a null test.

## Optional, if time allows

- **S-PrediXcan** (MetaXcan, GTEx v8 whole-blood models) on the ADHD GWAS gives predicted blood-expression associations for about 38,000 ADHD cases. The ME11 genes can then be tested as a set. This is a second platform-independent ADHD test.
- **Cabana-Domínguez summary statistics** (if shared): run `limma::cameraPR(stat_vector, list(ME11 = genes))` on their t-statistics, using the gene set from `ref$module_genes`.

# Changes to the revision scripts

These are copies of the scripts in `revision_analyses/`; the originals are unchanged.
Every change is marked `# PATCH:` in the code.

| File | Change | Why |
|---|---|---|
| `00_make_inputs.R` | **New.** Builds the five `inputs/` files from what the pipeline already saved (`ASD_wgcna_input.rds`, `ASD_metadata_clean.rds`, `ASD_network_final.rds`, the cached GEO series matrices and GPL10558 annotation). Also saves an all-probe pre-ComBat matrix. | Replaces README step 0, so the WGCNA script does not need editing or re-running. Stops unless it reproduces 238 samples, 115/123, the 28/69 and 87/54 split, 298 ME11 probes, 187 ME11 genes, 47,323 common probes and ME11 r = -0.26. |
| `00_config.R` | Adds `f_expr_precombat_full` and `max_block`; sets a CRAN mirror. | Used by 03 and 02; without a mirror `need()` cannot install missing packages. |
| `00_config.R` | `load_reference()` now defines each module's genes from **all** its probes (`module_sets`; ME11 = 187 genes, as in the paper). `module_genes` comes from this. | The one-probe-per-gene collapse keeps each gene's most highly expressed probe; for 87 of the 187 ME11 genes that probe is in another module, so ME11 shrank to 100 genes. The single-label `gene_labels` is still used for modulePreservation and the consensus overlap. |
| `04_sfari_brain_overlap.R`, `05a_make_magma_sets.R` | Use `module_sets` (all probes of each module). | Same reason. |
| `01_asd_replication.R` | Original cohorts re-scored with the same ME11 score as the new cohorts, instead of the paper's eigengene d values (-0.79, +0.27). | The meta-analysis compared two different measures (eigengene vs mean-z score); the eigengene sign is also arbitrary. |
| `01_asd_replication.R` | Effects reported for two ASD definitions: broad (autism + PDD-NOS + Asperger's) and strict (autism only). Adds a replication-only meta-analysis (new cohorts without the discovery cohorts). | PDD-NOS was excluded from GSE111175 in the paper. The replication-only estimate is the cleanest test. |
| `01_asd_replication.R` | `\btypical` instead of `typical` in the diagnosis parser; GEO downloads cached in `geo_cache/`; longer download timeout; skip cohorts with only one group. | "atypical autism" was read as a control; reruns no longer download again. |
| `01_asd_replication.R` | Mothers in GSE26415 ("healthy women who had children with ASD (asdMO)", "control (ctrlMO)") are excluded. Download retried 3 times with the reason printed; a hung download gives up after 15 min. Cohorts can be re-run one at a time (`Rscript 01_asd_replication.R GSE18123`) and earlier results are kept. | The mothers were counted as 21 extra ASD cases and 21 extra controls, which made the GSE26415 result invalid. GSE18123 and GSE25507 failed to download on the first run without saying why. |
| `01_asd_replication.R` | Diagnosis read only from GEO's parsed diagnosis/group field when one exists. | In GSE18123 every sample's description mentions autism, so all 82 controls were called ASD and the cohort was skipped. |
| `02_consensus_wgcna.R` | `maxBlockSize = max_block` (8000) instead of 20000. | One 15,000-gene block for two datasets needs about 6 GB RAM. |
| `03_cell_composition.R` | xCell uses all 47,323 probes (pre-ComBat), not only the 15,000 network probes. `moduleEigengenes()` given its required `colors` argument. | xCell scores about 10,800 signature genes; many are not among the most variable probes. The missing argument stopped the script. |
| `05b_run_magma.sh` | Checks that MAGMA starts and runs it through Rosetta on Apple-silicon Macs if needed. Reads the GWAS files with `gzip -dcf` instead of `zcat -f`. | The Mac binary may be Intel-only; macOS `zcat` cannot read `.gz` files. |
| `06_ppi_permutation.R` | `moduleEigengenes()` given its required `colors` argument; creates `string_cache/` and sets a longer download timeout. | The missing argument stopped the script; STRINGdb fails if its download folder does not exist. |


## New scripts (not in the original set)

| File | What it does |
|---|---|
| `07_consensus_modules_replication.R` | Tests every consensus module from 02 (and the original ME11) in all cohorts with the same score and labels as 01; FDR across modules. |
| `08_probe_quality.R` | Probe quality of ME11 vs other modules (intensity, secondary probes, annotation) and the ME11-ASD effect by probe subset. |

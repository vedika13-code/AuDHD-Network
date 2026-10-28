#!/usr/bin/env bash
# =====================================================================
# 05b_run_magma.sh
# Tests whether ME11 (and every other module) carries ASD or ADHD
# common-variant risk, using public GWAS summary statistics.
# This is the platform-independent ADHD test (reviewer point 1).
#
# Before running, download into magma/ (see README for links):
#   magma                      MAGMA v1.10 binary (CTG lab, VU Amsterdam)
#   g1000_eur.{bed,bim,fam}    1000 Genomes EUR reference (from the MAGMA site)
#   NCBI37.3.gene.loc          gene locations, build 37 (from the MAGMA site)
#   asd_gwas.txt               PGC ASD GWAS (Grove et al. 2019, iPSYCH-PGC)
#   adhd_gwas.txt              PGC ADHD GWAS (Demontis et al. 2023)
# Both GWAS are on hg19/GRCh37, matching NCBI37.3.
# Outputs: magma/asd_modules.gsa.out, magma/adhd_modules.gsa.out
# =====================================================================
set -euo pipefail
cd magma

# PATCH: on Apple-silicon Macs the MAGMA binary may be Intel-only; run it through Rosetta if needed
chmod +x magma
runs () { set +e; "$@" >/dev/null 2>&1; rc=$?; set -e; [ $rc -ne 126 ] && [ $rc -ne 127 ] && [ $rc -lt 128 ]; }  # 126/127 = cannot execute
if runs ./magma; then MAGMA=./magma
elif runs arch -x86_64 ./magma; then MAGMA="arch -x86_64 ./magma"
else echo "MAGMA does not start. On a Mac run: softwareupdate --install-rosetta --agree-to-license, and check that macOS did not block it (System Settings > Privacy & Security > Allow)."; exit 1; fi
echo "Using: $MAGMA"

# Sample sizes (cases + controls). Check against the readme of the file you download.
N_ASD=46350      # 18,381 cases + 27,969 controls (Grove 2019)
N_ADHD=225534    # 38,691 cases + 186,843 controls (Demontis 2023)

# Reduce a GWAS file to SNP and P columns, whatever order the columns are in
# PATCH: gzip -dcf instead of zcat -f (macOS zcat cannot read .gz); reads plain or gzipped files
prep () {  # $1 input, $2 output
  awk 'NR==1{for(i=1;i<=NF;i++){h=toupper($i); if(h=="SNP"||h=="RSID"||h=="MARKERNAME")s=i; if(h=="P"||h=="PVAL"||h=="P_VALUE")p=i}
            if(!s||!p){print "SNP or P column not found in header" > "/dev/stderr"; exit 1}
            print "SNP\tP"; next}
       {print $s"\t"$p}' <(gzip -dcf "$1") > "$2"
}
prep asd_gwas.txt  asd_snp_p.txt
prep adhd_gwas.txt adhd_snp_p.txt

# 1) map SNPs to genes (35 kb upstream, 10 kb downstream)
$MAGMA --annotate window=35,10 --snp-loc g1000_eur.bim --gene-loc NCBI37.3.gene.loc --out annot

# 2) gene-level association
$MAGMA --bfile g1000_eur --pval asd_snp_p.txt  N=$N_ASD  --gene-annot annot.genes.annot --out asd
$MAGMA --bfile g1000_eur --pval adhd_snp_p.txt N=$N_ADHD --gene-annot annot.genes.annot --out adhd

# 3) competitive gene-set analysis of the modules
#    a) against all genes in the genome
$MAGMA --gene-results asd.genes.raw  --set-annot module_sets.txt --out asd_modules
$MAGMA --gene-results adhd.genes.raw --set-annot module_sets.txt --out adhd_modules
#    b) conditioned on NETWORK_ALL, so modules are compared with other blood-expressed genes
#       (check the --model condition syntax in the MAGMA v1.10 manual if this line errors)
$MAGMA --gene-results asd.genes.raw  --set-annot module_sets.txt --model condition-hide=NETWORK_ALL --out asd_modules_cond
$MAGMA --gene-results adhd.genes.raw --set-annot module_sets.txt --model condition-hide=NETWORK_ALL --out adhd_modules_cond

echo "Done. Send back: magma/*_modules*.gsa.out"

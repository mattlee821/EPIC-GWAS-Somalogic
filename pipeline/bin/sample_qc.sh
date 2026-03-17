#!/bin/bash
set -euo pipefail

POSITIONAL_ARGS=()
bfile=""
study=""
mind=""
king=""
keep=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --bfile) bfile="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --mind) mind="$2"; shift 2 ;;
    --king) king="$2"; shift 2 ;;
    --keep) keep="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

# Detect bfile vs pfile
plink_flag="--bfile"
if [[ -f "${bfile}.pgen" ]]; then
  plink_flag="--pfile"
elif [[ ! -f "${bfile}.bed" ]]; then
  echo "ERROR: Neither ${bfile}.pgen nor ${bfile}.bed found."
  exit 1
fi

mkdir -p "$outdir"
log_file="${outdir}/sample_qc.log"
exec > >(tee -a "$log_file") 2>&1

echo "=========================================="
echo " Starting sample QC for $study"
echo "=========================================="

KEEP_ARG=""
if [[ -n "$keep" && -f "$keep" ]]; then
    echo "Restricting analysis to samples in: $keep"
    KEEP_ARG="--keep $keep"
fi

# Step 1: Mind filter
echo ">>> Step 1: Call rate filter (--mind $mind)"
plink2 $plink_flag "$bfile" \
  $KEEP_ARG \
  --mind "$mind" \
  --make-just-fam \
  --out "${outdir}/mind_pass"

# Step 2: Heterozygosity
echo ">>> Step 2: Heterozygosity check"
plink2 $plink_flag "$bfile" \
  --keep "${outdir}/mind_pass.fam" \
  --indep-pairwise 200 50 0.25 \
  --out "${outdir}/het_ld_temp"

plink2 $plink_flag "$bfile" \
  --keep "${outdir}/mind_pass.fam" \
  --extract "${outdir}/het_ld_temp.prune.in" \
  --het \
  --out "${outdir}/het"

# Use R for robust filtering (better than awk for variable columns)
echo "Filtering heterozygosity outliers (mean +/- 3SD)..."
Rscript - <<EOF
library(data.table)
het <- fread("${outdir}/het.het")
# Find F column dynamically
f_col <- grep("^F$", colnames(het), value=TRUE)
if(length(f_col) == 0) f_col <- colnames(het)[ncol(het)] # fallback to last col

f_vals <- het[[f_col]]
m <- mean(f_vals, na.rm=TRUE)
s <- sd(f_vals, na.rm=TRUE)
lower <- m - 3*s
upper <- m + 3*s

cat(sprintf("F-coeff mean: %.4f, SD: %.4f\n", m, s))
cat(sprintf("Keep range: [%.4f, %.4f]\n", lower, upper))

# Identify ID column (IID or #IID)
id_col <- grep("IID", colnames(het), value=TRUE)[1]
pass_ids <- het[get(f_col) >= lower & get(f_col) <= upper, ..id_col]
fwrite(pass_ids, "${outdir}/het_pass.txt", col.names=FALSE, sep="\t")
EOF

# Step 3: King cutoff
echo ">>> Step 3: Relatedness filter (--king-cutoff $king)"
plink2 $plink_flag "$bfile" \
  --keep "${outdir}/het_pass.txt" \
  --king-cutoff "$king" \
  --out "${outdir}/king_pass"

# Final clean sample list
if [[ -f "${outdir}/king_pass.king.cutoff.in.id" ]]; then
    cp "${outdir}/king_pass.king.cutoff.in.id" "${outdir}/qc_pass_samples.txt"
else
    # In case no relatives found to prune
    cp "${outdir}/het_pass.txt" "${outdir}/qc_pass_samples.txt"
fi

final_count=$(wc -l < "${outdir}/qc_pass_samples.txt")
echo "=========================================="
echo " QC Finished: $final_count samples passed all filters."
echo " Final list: ${outdir}/qc_pass_samples.txt"
echo "=========================================="

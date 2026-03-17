#!/bin/bash
set -euo pipefail

bfile=""
keep=""
extract=""
study=""
window=""
step=""
r2=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --bfile) bfile="$2"; shift 2 ;;
    --keep) keep="$2"; shift 2 ;;
    --extract) extract="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --window) window="$2"; shift 2 ;;
    --step) step="$2"; shift 2 ;;
    --r2) r2="$2"; shift 2 ;;
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
log_file="${outdir}/ld_pruning.log"
exec > >(tee -a "$log_file") 2>&1

echo "LD Pruning for $study using $plink_flag $bfile..."

plink2 $plink_flag "$bfile" \
  --keep "$keep" \
  --extract "$extract" \
  --indep-pairwise "${window}kb" "$step" "$r2" \
  --out "${outdir}/ld_pruned"

plink2 $plink_flag "$bfile" \
  --keep "$keep" \
  --extract "${outdir}/ld_pruned.prune.in" \
  --make-bed \
  --out "${outdir}/step1_input"

echo "Pruned step1 input files written to ${outdir}/step1_input.*"

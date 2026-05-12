#!/bin/bash
set -euo pipefail

pfile=""
keep=""
extract=""
study=""
window=""
step=""
r2=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --pfile) pfile="$2"; shift 2 ;;
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

if [[ -z "$pfile" || -z "$keep" || -z "$extract" || -z "$study" || -z "$window" || -z "$step" || -z "$r2" || -z "$outdir" ]]; then
  echo "Usage: $0 --pfile <prefix> --keep <file> --extract <file> --study <val> --window <kb> --step <val> --r2 <val> --outdir <dir>" >&2
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"
log_file="${outdir}/ld_pruning.log"
exec > >(tee -a "$log_file") 2>&1

echo "LD Pruning for $study using PLINK2 prefix $pfile..."

plink2 --pfile "$pfile" \
  --keep "$keep" \
  --extract "$extract" \
  --indep-pairwise "${window}kb" "$step" "$r2" \
  --out "${outdir}/ld_pruned"

plink2 --pfile "$pfile" \
  --keep "$keep" \
  --extract "${outdir}/ld_pruned.prune.in" \
  --make-pgen \
  --out "${outdir}/step1_input"

echo "Pruned step1 input files written to ${outdir}/step1_input.*"

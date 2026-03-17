#!/bin/bash
set -euo pipefail

bfile=""
keep=""
study=""
maf=""
hwe=""
geno=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --bfile) bfile="$2"; shift 2 ;;
    --keep) keep="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --maf) maf="$2"; shift 2 ;;
    --hwe) hwe="$2"; shift 2 ;;
    --geno) geno="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$bfile" || -z "$keep" || -z "$study" || -z "$maf" || -z "$hwe" || -z "$geno" || -z "$outdir" ]]; then
  echo "Usage: $0 --bfile <val> --keep <val> --study <val> --maf <val> --hwe <val> --geno <val> --outdir <val>"
  exit 1
fi

# Detect bfile vs pfile
plink_flag="--bfile"
if [[ -f "${bfile}.pgen" ]]; then
  plink_flag="--pfile"
elif [[ ! -f "${bfile}.bed" ]]; then
  echo "ERROR: Neither ${bfile}.pgen nor ${bfile}.bed found."
  exit 1
fi

mkdir -p "$outdir"
log_file="${outdir}/variant_qc.log"
exec > >(tee -a "$log_file") 2>&1

echo "Starting variant QC for $study using $plink_flag $bfile..."

plink2 $plink_flag "$bfile" \
  --keep "$keep" \
  --maf "$maf" \
  --hwe "$hwe" \
  --geno "$geno" \
  --write-snplist \
  --out "${outdir}/variant_qc_pass"

echo "Variant QC complete. Variant list written to ${outdir}/variant_qc_pass.snplist"

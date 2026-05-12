#!/bin/bash
set -euo pipefail

pfile=""
keep=""
study=""
maf=""
hwe=""
geno=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --pfile) pfile="$2"; shift 2 ;;
    --keep) keep="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --maf) maf="$2"; shift 2 ;;
    --hwe) hwe="$2"; shift 2 ;;
    --geno) geno="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$pfile" || -z "$keep" || -z "$study" || -z "$maf" || -z "$hwe" || -z "$geno" || -z "$outdir" ]]; then
  echo "Usage: $0 --pfile <prefix> --keep <val> --study <val> --maf <val> --hwe <val> --geno <val> --outdir <val>"
  exit 1
fi

for ext in pgen pvar psam; do
  if [[ ! -f "${pfile}.${ext}" ]]; then
    echo "ERROR: Missing PLINK2 input file: ${pfile}.${ext}" >&2
    exit 1
  fi
done

mkdir -p "$outdir"
log_file="${outdir}/variant_qc.log"
exec > >(tee -a "$log_file") 2>&1

echo "Starting variant QC for $study using PLINK2 prefix $pfile..."

plink2 --pfile "$pfile" \
  --keep "$keep" \
  --maf "$maf" \
  --hwe "$hwe" \
  --geno "$geno" \
  --write-snplist \
  --out "${outdir}/variant_qc_pass"

echo "Variant QC complete. Variant list written to ${outdir}/variant_qc_pass.snplist"

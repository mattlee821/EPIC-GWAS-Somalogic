#!/bin/bash
set -euo pipefail

step1_bfile=""
pheno_file=""
cov_file=""
study=""
group=""
bsize=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --step1_bfile) step1_bfile="$2"; shift 2 ;;
    --pheno_file) pheno_file="$2"; shift 2 ;;
    --cov_file) cov_file="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --bsize) bsize="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"
mkdir -p "${outdir}/tmp"

log_file="${outdir}/step1.log"

{
  echo "COVARIATE FILE PREVIEW (first 2 lines):"
  head -n 2 "${cov_file}" 2>/dev/null || echo "WARNING: failed to read covariate file '${cov_file}'"
  echo ""
} > "$log_file"

if ! regenie \
  --step 1 \
  --bed "${step1_bfile}" \
  --phenoFile "${pheno_file}" \
  --covarFile "${cov_file}" \
  --bsize "${bsize}" \
  --lowmem \
  --lowmem-prefix "${outdir}/tmp/rg" \
  --threads 8 \
  --gz \
  --force-step1 \
  --out "${outdir}/step1" >> "$log_file" 2>&1; then
    echo "ERROR: REGENIE Step 1 failed. Log follows:"
    cat "$log_file"
    exit 1
fi

if [[ -f "${outdir}/step1_pred.list" ]]; then
  mv "${outdir}/step1_pred.list" "${outdir}/pred.list"
fi

for f in "${outdir}"/step1_*.loco.gz; do
  [[ -e "$f" ]] || continue
  base="${f##*/step1_}"
  mv "$f" "${outdir}/loco_${base}"
done

# Update pred.list to match renamed LOCO files
if [[ -f "${outdir}/pred.list" ]]; then
  sed -i 's/step1_/loco_/g' "${outdir}/pred.list"
fi

echo "REGENIE Step 1 for $study ($group) complete."

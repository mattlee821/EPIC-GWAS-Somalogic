#!/bin/bash
set -euo pipefail

bgen_file=""
bgen_dir=""
pgen_prefix=""
sample_file=""
pheno_file=""
cov_file=""
pred_list=""
study=""
group=""
chromosomes=""
bsize=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --bgen_file) bgen_file="$2"; shift 2 ;;
    --bgen_dir) bgen_dir="$2"; shift 2 ;;
    --pfile) pgen_prefix="$2"; shift 2 ;;
    --sample_file) sample_file="$2"; shift 2 ;;
    --pheno_file) pheno_file="$2"; shift 2 ;;
    --cov_file) cov_file="$2"; shift 2 ;;
    --pred_list) pred_list="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --chr|--chromosomes) chromosomes="$2"; shift 2 ;;
    --bsize) bsize="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"

if [[ -z "$chromosomes" ]]; then
    echo "ERROR: At least one chromosome must be provided via --chromosomes"
    exit 1
fi

# Default to Nextflow-provided CPU count when available (e.g., NXF_TASK_CPUS).
# Fall back to SLURM_CPUS_PER_TASK, then 1 for manual runs.
threads="${NXF_TASK_CPUS:-${SLURM_CPUS_PER_TASK:-1}}"
export OMP_NUM_THREADS="${threads}"
export OPENBLAS_NUM_THREADS="${threads}"
export MKL_NUM_THREADS="${threads}"
export NUMEXPR_MAX_THREADS="${threads}"

if [[ -n "$pred_list" ]]; then
    localized_pred_list="${outdir}/pred.list"
    python3 - "$pred_list" "$localized_pred_list" <<'PY'
import os
import re
import sys

src_path, dst_path = sys.argv[1:]
tmp_path = dst_path + ".tmp"
pattern = re.compile(r'\S+\.loco\.gz')

def rewrite_token(match: re.Match[str]) -> str:
    token = match.group(0)
    base = os.path.basename(token)
    if base.endswith(".loco.gz") and os.path.exists(base):
        return f"./{base}"
    return token

with open(src_path, "r", encoding="utf-8") as src, open(tmp_path, "w", encoding="utf-8") as dst:
    for raw_line in src:
        dst.write(pattern.sub(rewrite_token, raw_line))

os.replace(tmp_path, dst_path)
PY
    pred_list="${localized_pred_list}"
fi

if [[ -n "$pgen_prefix" ]]; then
    # REGENIE 4.x check: It requires #FID as 1st col and IID as 2nd col.
    ORIG_PSAM="${pgen_prefix}.psam"
    header=$(head -n 1 "$ORIG_PSAM")

    if [[ "$header" != *"#FID"* ]]; then
        echo ">>> Fixing PSAM header for REGENIE compatibility..."

        # Create a new prefix in the output directory
        FIXED_PREFIX="${outdir}/fixed_gen"

        # 1. Create fixed PSAM with strict #FID IID format
        # If the input starts with #IID, we strip the # from IID and prepend #FID
        awk 'BEGIN {OFS="\t"} 
             NR==1 {
                sub(/^#/,"",$1); 
                print "#FID", "IID", "SEX", "PHENO"
             } 
             NR>1 {
                print "0", $1, $2, $3
             }' "$ORIG_PSAM" > "${FIXED_PREFIX}.psam"

        # 2. Symlink the large data files to this new prefix
        ln -sf "$(readlink -f ${pgen_prefix}.pgen)" "${FIXED_PREFIX}.pgen"
        ln -sf "$(readlink -f ${pgen_prefix}.pvar)" "${FIXED_PREFIX}.pvar"

        pgen_prefix="${FIXED_PREFIX}"
    fi
fi

IFS=',' read -r -a chr_array <<< "$chromosomes"

for chr in "${chr_array[@]}"; do
    chr="${chr//[[:space:]]/}"
    [[ -n "$chr" ]] || continue

    log_file="${outdir}/step2_chr${chr}.log"
    cmd=(regenie --step 2 --threads "${threads}")

    if [[ -n "$bgen_dir" ]] && [[ -f "${bgen_dir}/chr${chr}.bgen" ]]; then
        cmd+=(--bgen "${bgen_dir}/chr${chr}.bgen" --sample "${sample_file}")
    elif [[ -n "$bgen_file" ]]; then
        cmd+=(--bgen "${bgen_file}" --sample "${sample_file}")
    elif [[ -n "$pgen_prefix" ]]; then
        cmd+=(--pgen "${pgen_prefix}")
    else
        echo "ERROR: No genotype input available for chromosome ${chr}" | tee "$log_file"
        exit 1
    fi

    cmd+=(
        --phenoFile "${pheno_file}"
        --covarFile "${cov_file}"
        --pred "${pred_list}"
        --chr "${chr}"
        --bsize "${bsize}"
        --qt
        --gz
        --out "${outdir}/full_chr${chr}"
    )

    printf 'Running: %q ' "${cmd[@]}" | tee "$log_file"
    printf '\n' | tee -a "$log_file"

    if ! "${cmd[@]}" >> "$log_file" 2>&1; then
        echo "ERROR: REGENIE Step 2 failed. Log follows:"
        cat "$log_file"
        exit 1
    fi
done

echo "REGENIE Step 2 for $study ($group) complete across chromosomes: $chromosomes"

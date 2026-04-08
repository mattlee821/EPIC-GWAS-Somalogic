#!/bin/bash
set -euo pipefail

bgen_file=""
pgen_prefix=""
sample_file=""
pheno_file=""
cov_file=""
pred_list=""
study=""
group=""
chr=""
bsize=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --bgen_file) bgen_file="$2"; shift 2 ;;
    --pfile) pgen_prefix="$2"; shift 2 ;;
    --sample_file) sample_file="$2"; shift 2 ;;
    --pheno_file) pheno_file="$2"; shift 2 ;;
    --cov_file) cov_file="$2"; shift 2 ;;
    --pred_list) pred_list="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --chr) chr="$2"; shift 2 ;;
    --bsize) bsize="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"
log_file="${outdir}/step2_chr${chr}.log"

# Default to Nextflow-provided CPU count when available (e.g., NXF_TASK_CPUS).
# Fall back to SLURM_CPUS_PER_TASK, then 2 for manual runs.
threads="${NXF_TASK_CPUS:-${SLURM_CPUS_PER_TASK:-2}}"

# Build REGENIE command
CMD="regenie --step 2 --threads ${threads}"

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
        
        CMD+=" --pgen ${FIXED_PREFIX}"
    else
        CMD+=" --pgen ${pgen_prefix}"
    fi
elif [[ -n "$bgen_file" ]]; then
    CMD+=" --bgen ${bgen_file} --sample ${sample_file}"
fi

CMD+=" --phenoFile ${pheno_file}"
CMD+=" --covarFile ${cov_file}"
CMD+=" --pred ${pred_list}"
CMD+=" --chr ${chr}"
CMD+=" --bsize ${bsize}"
CMD+=" --qt --gz"
CMD+=" --out ${outdir}/full_chr${chr}"

echo "Running: $CMD"
if ! eval "$CMD" > "$log_file" 2>&1; then
    echo "ERROR: REGENIE Step 2 failed. Log follows:"
    cat "$log_file"
    exit 1
fi

echo "REGENIE Step 2 for $study ($group) chr$chr complete."

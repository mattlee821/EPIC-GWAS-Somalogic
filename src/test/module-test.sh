#!/bin/bash
# Script: src/test/module-test.sh
# Purpose: Manually execute individual pipeline modules for debugging.
# This script handles environment switching automatically where possible.

set -euo pipefail

# Helper to activate conda within the script
# We try to source the standard conda profile.
# If this fails, the user must activate the environment manually before running.
CONDA_PATH=$(conda info --base 2>/dev/null || echo "")
if [[ -f "${CONDA_PATH}/etc/profile.d/conda.sh" ]]; then
    source "${CONDA_PATH}/etc/profile.d/conda.sh"
else
    echo "Warning: Could not find conda.sh. Automatic environment switching may fail."
fi

# Prevent R_HOME conflicts
unset R_HOME

# --- Global Test Variables ---
PROJECT_ROOT=$(pwd)
PIPELINE_BIN="${PROJECT_ROOT}/pipeline/bin"
DATA_DIR="${PROJECT_ROOT}/data"
TEST_OUT="${PROJECT_ROOT}/test_run"
mkdir -p "${TEST_OUT}"

PHENO_FILE="${DATA_DIR}/phenotype/phenofile-test.txt"
SAMPLESHEET="${PROJECT_ROOT}/pipeline/samplesheet.csv"
STUDY="Neuro_01"
# For Neuro_01, we know it is pfile (PLINK2)
GENETIC_PREFIX="${DATA_DIR}/genetics/Neuro_01/EPIC_GSA_imputed_QC_plink"

usage() {
    echo "Usage: $0 [module_name]"
    echo "Available modules:"
    echo "  validate    - Run input validation (nf_master)"
    echo "  pheno_prep  - Run phenotype preparation (r_gwas)"
    echo "  cov_prep    - Run covariate preparation (r_gwas)"
    echo "  sample_qc   - Run Sample QC (plink2)"
    echo "  variant_qc  - Run Variant QC (plink2)"
    echo "  ld_prun     - Run LD Pruning (plink2)"
    exit 1
}

if [[ $# -lt 1 ]]; then usage; fi

MODULE=$1

case $MODULE in
    validate)
        echo ">>> Testing VALIDATE_INPUTS (using nf_master)..."
        conda activate nf_master 2>/dev/null || echo "Info: nf_master not active, proceeding with current env"
        python3 "${PIPELINE_BIN}/validate_inputs.py" \
            --samplesheet "${SAMPLESHEET}" \
            --phenotype_file "${PHENO_FILE}" \
            --covariate_file "${PHENO_FILE}" \
            --include_studies "${STUDY}" \
            --outdir "${TEST_OUT}/validate"
        ;;

    pheno_prep)
        echo ">>> Testing PREPARE_PHENOTYPES (using r_gwas)..."
        conda activate r_gwas 2>/dev/null || echo "Info: r_gwas not active/found, using current R"
        Rscript "${PIPELINE_BIN}/prepare_phenotypes.R" \
            --phenotype_file "${PHENO_FILE}" \
            --sample_file "${GENETIC_PREFIX}.psam" \
            --study_id "${STUDY}" \
            --group "all" \
            --covariate_file "${PHENO_FILE}" \
            --include_proteins "10620-21,18878-15" \
            --chunk_size 10 \
            --outdir "${TEST_OUT}/pheno_prep"
        ;;

    cov_prep)
        echo ">>> Testing PREPARE_COVARIATES (using r_gwas)..."
        conda activate r_gwas 2>/dev/null || echo "Info: r_gwas not active/found, using current R"
        Rscript "${PIPELINE_BIN}/prepare_covariates.R" \
            --covariate_file "${PHENO_FILE}" \
            --sample_file "${GENETIC_PREFIX}.psam" \
            --study_id "${STUDY}" \
            --group "all" \
            --include_covariates "Sex,Age_Recr" \
            --outdir "${TEST_OUT}/cov_prep"
        ;;

    sample_qc)
        echo ">>> Testing SAMPLE_QC (using plink2)..."
        conda activate plink2 2>/dev/null || { echo "ERROR: plink2 environment not found."; exit 1; }
        
        # 1. Create a restricted sample list from the matched pheno prep
        PHENO_PREP="${TEST_OUT}/pheno_prep/${STUDY}_all_full.pheno"
        RESTRICTED_LIST="${TEST_OUT}/sample_qc/matched_samples.txt"
        mkdir -p "${TEST_OUT}/sample_qc"
        
        if [ ! -f "$PHENO_PREP" ]; then
            echo "Error: Run pheno_prep first to define matched samples."
            exit 1
        fi
        
        # Extract only IID (column 2) from the pheno file for PLINK2 keep list
        awk 'NR>1 {print $2}' "$PHENO_PREP" > "$RESTRICTED_LIST"
        echo "Restricting QC to $(wc -l < "$RESTRICTED_LIST") matched samples."

        bash "${PIPELINE_BIN}/sample_qc.sh" \
            --bfile "${GENETIC_PREFIX}" \
            --keep "$RESTRICTED_LIST" \
            --study "${STUDY}" \
            --mind 0.05 \
            --king 0.0884 \
            --outdir "${TEST_OUT}/sample_qc"
        ;;

    variant_qc)
        echo ">>> Testing VARIANT_QC (using plink2)..."
        conda activate plink2 2>/dev/null || { echo "ERROR: plink2 environment not found."; exit 1; }
        SAMPLE_LIST="${TEST_OUT}/sample_qc/${STUDY}_qc_pass_samples.txt"
        if [ ! -f "$SAMPLE_LIST" ]; then echo "Error: Run sample_qc first"; exit 1; fi
        
        bash "${PIPELINE_BIN}/variant_qc.sh" \
            --bfile "${GENETIC_PREFIX}" \
            --keep "$SAMPLE_LIST" \
            --study "${STUDY}" \
            --maf 0.01 \
            --hwe 1e-15 \
            --geno 0.05 \
            --outdir "${TEST_OUT}/variant_qc"
        ;;

    ld_prun)
        echo ">>> Testing LD_PRUNING (using plink2)..."
        conda activate plink2 2>/dev/null || { echo "ERROR: plink2 environment not found."; exit 1; }
        SAMPLE_LIST="${TEST_OUT}/sample_qc/${STUDY}_qc_pass_samples.txt"
        VARIANT_LIST="${TEST_OUT}/variant_qc/${STUDY}_variant_qc_pass.snplist"
        
        bash "${PIPELINE_BIN}/ld_pruning.sh" \
            --bfile "${GENETIC_PREFIX}" \
            --keep "$SAMPLE_LIST" \
            --extract "$VARIANT_LIST" \
            --study "${STUDY}" \
            --window 1000 \
            --step 1 \
            --r2 0.5 \
            --outdir "${TEST_OUT}/ld_prun"
        ;;

    regenie_step1)
        echo ">>> Testing REGENIE Step 1 (using regenie env)..."
        conda activate regenie 2>/dev/null || { echo "ERROR: regenie environment not found."; exit 1; }
        
        PHENO="${TEST_OUT}/pheno_prep/${STUDY}_all_full.pheno"
        COVAR="${TEST_OUT}/cov_prep/${STUDY}_all.cov"
        PRUNED_BED_PREFIX="${TEST_OUT}/ld_prun/${STUDY}_step1_input"
        
        mkdir -p "${TEST_OUT}/regenie_step1"
        
        regenie --step 1 \
            --bed "$PRUNED_BED_PREFIX" \
            --phenoFile "$PHENO" \
            --covarFile "$COVAR" \
            --bsize 1000 \
            --lowmem \
            --gz \
            --threads 8 \
            --force-step1 \
            --out "${TEST_OUT}/regenie_step1/${STUDY}_step1"
        ;;

    regenie_step2)
        echo ">>> Testing REGENIE Step 2 (using regenie env)..."
        conda activate regenie 2>/dev/null || { echo "ERROR: regenie environment not found."; exit 1; }
        
        # Neuro_01 is pfile based
        GEN_PREFIX="${DATA_DIR}/genetics/${STUDY}/EPIC_GSA_imputed_QC_plink"
        PHENO="${TEST_OUT}/pheno_prep/${STUDY}_all_chunk_001.pheno"
        COVAR="${TEST_OUT}/cov_prep/${STUDY}_all.cov"
        PRED_LIST="${TEST_OUT}/regenie_step1/${STUDY}_step1_pred.list"
        
        # We still test one chromosome (21) for speed
        CHR="21"
        
        mkdir -p "${TEST_OUT}/regenie_step2"
        
        bash "${PIPELINE_BIN}/run_regenie_step2.sh" \
            --pfile "$GEN_PREFIX" \
            --pheno_file "$PHENO" \
            --cov_file "$COVAR" \
            --pred_list "$PRED_LIST" \
            --study "$STUDY" \
            --group "all" \
            --chunk_id "chunk_001" \
            --chr "$CHR" \
            --bsize 400 \
            --outdir "${TEST_OUT}/regenie_step2"
        ;;

    *)
        usage
        ;;
esac

echo "======================================"
echo "Test for $MODULE finished."
echo "Check outputs in: ${TEST_OUT}/${MODULE}"
echo "======================================"

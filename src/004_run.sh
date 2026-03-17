#!/bin/bash
# Script: src/004_run.sh
# Purpose: Wrapper to quickly launch the Nextflow GWAS pipeline with runtime overrides.

set -euo pipefail

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

# Initialize conda for the bash script subshell
source "$(conda info --base)/etc/profile.d/conda.sh" || true
conda activate nf_master
unset R_HOME

# --- Default Variables ---
PIPELINE_DIR="pipeline"
PARAMS_FILE="params.yaml"
PROFILE="slurm"
RESULT_DIR="results"

# Runtime overrides
INCLUDE_STUDIES=""
INCLUDE_PROTEINS=""
PHENO_FILE=""
COVAR_FILE=""
COVARIATES=""
META_GROUP=""
META_STUDY=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    --out|-o|--outdir)
      RESULT_DIR="$2"
      shift 2
      ;;
    --studies)
      INCLUDE_STUDIES="$2"
      shift 2
      ;;
    --proteins)
      INCLUDE_PROTEINS="$2"
      shift 2
      ;;
    --pheno)
      PHENO_FILE="$2"
      shift 2
      ;;
    --covar)
      COVAR_FILE="$2"
      shift 2
      ;;
    --covariates)
      COVARIATES="$2"
      shift 2
      ;;
    --meta-group)
      META_GROUP="$2"
      shift 2
      ;;
    --meta-study)
      META_STUDY="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --profile     Nextflow profile (default: slurm; alternatives: test, local)"
      echo "  --out, -o     Output directory for results (default: ./results)"
      echo "  --studies     Comma-separated list of studies to run (override params)"
      echo "  --proteins    Comma-separated list of proteins to run (override params)"
      echo "  --pheno       Path to phenotype file (overrides params.yaml)"
      echo "  --covar       Path to covariate file (overrides params.yaml)"
      echo "  --covariates  Comma-separated list of covariates (e.g. 'Age,Sex,PC1')"
      echo "  --meta-group  true/false: run within-study cases+controls meta"
      echo "  --meta-study  Comma-separated list of analyses for study-level meta (cases,controls,combined,meta). Default: false"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Ensure we're in the right place relative to the pipeline
if [[ ! -d "$PIPELINE_DIR" ]]; then
    echo "ERROR: Pipeline directory '$PIPELINE_DIR' not found."
    echo "       Please run this script from the GWAS project root."
    exit 1
fi

# Resolve RESULT_DIR to absolute path
ABS_RESULT_DIR="$(cd "$(dirname "$RESULT_DIR")" 2>/dev/null && pwd)/$(basename "$RESULT_DIR")" \
  || ABS_RESULT_DIR="$(pwd)/${RESULT_DIR}"
mkdir -p "$ABS_RESULT_DIR"

echo "=========================================="
echo " Launching Nextflow GWAS Pipeline"
echo "=========================================="
echo "Profile:  $PROFILE"
echo "Studies:  ${INCLUDE_STUDIES:-ALL}"
echo "Proteins: ${INCLUDE_PROTEINS:-ALL}"
echo "Output:   $ABS_RESULT_DIR"
echo "Meta-group: ${META_GROUP:-false}"
echo "Meta-study: ${META_STUDY:-false}"
echo "=========================================="

# Build base nextflow command
NF_CMD="nextflow run main.nf -profile $PROFILE -resume -params-file $PARAMS_FILE"
NF_CMD+=" --outdir \"$ABS_RESULT_DIR\""

# Append overrides if provided
if [[ -n "$INCLUDE_STUDIES" ]]; then
    NF_CMD+=" --include_studies \"$INCLUDE_STUDIES\""
fi

if [[ -n "$INCLUDE_PROTEINS" ]]; then
    # If it's a file, resolve to absolute path
    if [[ -f "$INCLUDE_PROTEINS" ]]; then
        ABS_PROTS="$(cd "$(dirname "$INCLUDE_PROTEINS")" && pwd)/$(basename "$INCLUDE_PROTEINS")"
        NF_CMD+=" --include_proteins \"$ABS_PROTS\""
    else
        NF_CMD+=" --include_proteins \"$INCLUDE_PROTEINS\""
    fi
fi

if [[ -n "${PHENO_FILE:-}" ]]; then
    ABS_PHENO="$(cd "$(dirname "$PHENO_FILE")" && pwd)/$(basename "$PHENO_FILE")"
    NF_CMD+=" --phenotype_file \"$ABS_PHENO\""
fi

if [[ -n "${COVAR_FILE:-}" ]]; then
    ABS_COVAR="$(cd "$(dirname "$COVAR_FILE")" && pwd)/$(basename "$COVAR_FILE")"
    NF_CMD+=" --covariate_file \"$ABS_COVAR\""
fi

if [[ -n "${COVARIATES:-}" ]]; then
    NF_CMD+=" --covariates \"$COVARIATES\""
fi

if [[ -n "${META_GROUP:-}" ]]; then
    NF_CMD+=" --meta_group \"$META_GROUP\""
fi

if [[ -n "${META_STUDY:-}" ]]; then
    NF_CMD+=" --meta_study \"$META_STUDY\""
fi

echo "Executing:"
echo "cd $PIPELINE_DIR && $NF_CMD"

# Move to pipeline directory and execute
cd "$PIPELINE_DIR"
eval $NF_CMD

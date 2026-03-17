#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=logs/test.log
#SBATCH --error=logs/test.err
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_p

set -euo pipefail

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

# Resolve paths
PROJ_ROOT="${GWAS_PROJECT_ROOT:-$(pwd)}"
ANALYSIS_OUT="${GWAS_ANALYSIS_DIR:-${PROJ_ROOT}/analysis}/test/"

# Standardize working directory to project root
cd "$PROJ_ROOT"

start_time=$(date +%s)

# Load conda and activate environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nf_master
unset R_HOME

bash src/004_run.sh \
  --profile slurm \
  --out "$ANALYSIS_OUT" \
  --studies "Neuro_01,Brea_01_Erneg" \
  --proteins "10620-21,18878-15,4721-54,8323-163" \
  --pheno "data/phenotype/phenofile-test.txt" \
  --covariates "Sex,Age_Recr" \
  --meta-group true \
  --meta-study "cases,controls,combined,meta"

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "Total runtime: $(date -u -d @$elapsed +%H:%M:%S)"

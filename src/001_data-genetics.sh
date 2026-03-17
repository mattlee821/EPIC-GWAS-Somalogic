#!/bin/bash

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

# Environment Variables
# We strictly require these to be set in .env or the environment
DIR_OUT="${GWAS_GENETICS_DIR:?Error: GWAS_GENETICS_DIR not set in .env}"
EPIC4ND_DIR="${GWAS_EPIC4ND_ARCHIVE:?Error: GWAS_EPIC4ND_ARCHIVE not set in .env}"

mkdir -p "$DIR_OUT"

echo "Starting Manual Documented Migration..."

# NEURO
echo "Processing: Neuro_01"
mkdir -p "$DIR_OUT/Neuro_01/"
rsync -rvP "$EPIC4ND_DIR/"*.{pgen,psam,pvar} "$DIR_OUT/Neuro_01/"

# BREA
echo "Processing: Brea_01_Erneg"
mkdir -p "$DIR_OUT/Brea_01_Erneg/"
rsync -rvP "$EPIC4ND_DIR/"*.{pgen,psam,pvar} "$DIR_OUT/Brea_01_Erneg/"

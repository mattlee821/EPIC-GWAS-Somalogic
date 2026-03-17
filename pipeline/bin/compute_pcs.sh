#!/bin/bash
set -euo pipefail

step1_bfile=""
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --bfile) step1_bfile="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

if [[ -z "$step1_bfile" || -z "$outdir" ]]; then
  echo "Usage: $0 --bfile <step1_input> --outdir <outdir>"
  exit 1
fi

mkdir -p "$outdir"

plink2 --bfile "$step1_bfile" \
  --pca 10 \
  --out "${outdir}/pca"

if [[ ! -f "${outdir}/pca.eigenvec" ]]; then
  echo "ERROR: PCA failed to produce pca.eigenvec" >&2
  exit 1
fi

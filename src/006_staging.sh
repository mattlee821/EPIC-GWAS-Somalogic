#!/bin/bash
#SBATCH --job-name=gwas-staging
#SBATCH --output=src/logs/staging.log
#SBATCH --error=src/logs/staging.err
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

# Script: src/006_staging.sh
# Purpose: Stage across-study meta-analysis files for transfer to collaborators.

set -euo pipefail

SCRIPT_SOURCE="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(dirname -- "$SCRIPT_SOURCE")" && pwd)"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  DEFAULT_PROJ_ROOT="${SLURM_SUBMIT_DIR}"
else
  DEFAULT_PROJ_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi

source_env_file() {
  local env_file="$1"
  if [[ -f "$env_file" ]]; then
    set -a
    source "$env_file"
    set +a
  fi
}

resolve_path() {
  local path="$1"
  case "$path" in
    /*) printf '%s\n' "$path" ;;
    *) printf '%s\n' "${PROJ_ROOT%/}/$path" ;;
  esac
}

usage() {
  cat <<'EOF'
Usage:
  bash src/006_staging.sh [options]
  sbatch src/006_staging.sh [options]

Finds across-study meta-analysis files with this pipeline layout:
  <outdir>/.../meta/<feature>/013_QC/meta.tsv.gz

Stages each file with the feature directory name:
  <transfer-dir>/<feature>.tsv.gz

Options:
  --out, --outdir PATH       Analysis output root to search
                             (default: $GWAS_ANALYSIS_DIR or <project>/analysis)
  --transfer, --transfer-dir PATH
                             Staging directory (default: <project>/transfer)
  --dry-run                  Show what would be staged without copying files
  -h, --help                 Show this help
EOF
}

source_env_file "${DEFAULT_PROJ_ROOT}/.env"
PROJ_ROOT="${GWAS_PROJECT_ROOT:-$DEFAULT_PROJ_ROOT}"
mkdir -p "${PROJ_ROOT}/src/logs"

RESULT_DIR="${GWAS_ANALYSIS_DIR:-${PROJ_ROOT}/analysis}"
TRANSFER_DIR="${PROJ_ROOT}/transfer"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out|--outdir)
      RESULT_DIR="$2"
      shift 2
      ;;
    --transfer|--transfer-dir)
      TRANSFER_DIR="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

ABS_RESULT_DIR="$(resolve_path "$RESULT_DIR")"
ABS_TRANSFER_DIR="$(resolve_path "$TRANSFER_DIR")"

if [[ ! -d "$ABS_RESULT_DIR" ]]; then
  echo "ERROR: Analysis output directory not found: $ABS_RESULT_DIR" >&2
  exit 1
fi

if ! command -v rsync >/dev/null 2>&1; then
  echo "ERROR: rsync is required but was not found in PATH." >&2
  exit 1
fi

declare -a FILES=()
while IFS= read -r -d '' file; do
  FILES+=("$file")
done < <(
  find "$ABS_RESULT_DIR" \
    -path "$ABS_TRANSFER_DIR" -prune -o \
    \( -path '*/meta/*/013_QC/meta.tsv.gz' -o -path '*/meta/*/*/013_QC/meta.tsv.gz' \) \
    \( -type f -o -type l \) \
    -print0
)

if [[ "${#FILES[@]}" -eq 0 ]]; then
  echo "ERROR: No study-level meta-analysis files found under: $ABS_RESULT_DIR" >&2
  echo "Expected files matching: */meta/<feature>/013_QC/meta.tsv.gz" >&2
  exit 1
fi

declare -a DEST_KEYS=()
declare -a DEST_SRCS=()
declare -a PLAN_SRCS=()
declare -a PLAN_DESTS=()
declare -a PLAN_FEATURES=()
declare -a PLAN_GROUPS=()
BROKEN_LINK_COUNT=0

for src in "${FILES[@]}"; do
  qc_dir="$(dirname -- "$src")"
  parent_dir="$(dirname -- "$qc_dir")"
  grandparent_dir="$(dirname -- "$parent_dir")"
  greatgrandparent_dir="$(dirname -- "$grandparent_dir")"

  if [[ "$(basename -- "$qc_dir")" != "013_QC" ]]; then
    echo "WARNING: Skipping unexpected path: $src" >&2
    continue
  fi

  if [[ -L "$src" && ! -e "$src" ]]; then
    BROKEN_LINK_COUNT=$((BROKEN_LINK_COUNT + 1))
    echo "ERROR: Published meta output is a broken symlink: $src" >&2
    echo "Target: $(readlink -- "$src" || true)" >&2
    continue
  fi

  if [[ "$(basename -- "$grandparent_dir")" == "meta" ]]; then
    group="meta"
    feature="$(basename -- "$parent_dir")"
  elif [[ "$(basename -- "$greatgrandparent_dir")" == "meta" ]]; then
    group="$(basename -- "$parent_dir")"
    feature="$(basename -- "$grandparent_dir")"
  else
    echo "WARNING: Skipping unexpected path: $src" >&2
    continue
  fi

  rel_dest="${feature}.tsv.gz"
  dest="${ABS_TRANSFER_DIR}/${rel_dest}"

  for idx in "${!DEST_KEYS[@]}"; do
    if [[ "${DEST_KEYS[$idx]}" == "$rel_dest" ]]; then
      echo "ERROR: Multiple source files would stage to the same destination: $rel_dest" >&2
      echo "First:  ${DEST_SRCS[$idx]}" >&2
      echo "Second: $src" >&2
      echo "Use separate transfer directories or remove duplicate batch outputs before staging." >&2
      exit 1
    fi
  done
  DEST_KEYS+=("$rel_dest")
  DEST_SRCS+=("$src")
  PLAN_SRCS+=("$src")
  PLAN_DESTS+=("$dest")
  PLAN_FEATURES+=("$feature")
  PLAN_GROUPS+=("$group")
done

if [[ "$BROKEN_LINK_COUNT" -gt 0 ]]; then
  echo "ERROR: Found $BROKEN_LINK_COUNT broken published meta output symlink(s)." >&2
  echo "These files were published by Nextflow as symlinks, but their targets are not accessible from this node." >&2
  echo "Run staging on a node where the Nextflow work directory is mounted, or rerun/republish the pipeline with real copied outputs." >&2
  exit 1
fi

if [[ "${#PLAN_SRCS[@]}" -eq 0 ]]; then
  echo "ERROR: No valid study-level meta-analysis files found under: $ABS_RESULT_DIR" >&2
  exit 1
fi

if [[ "$DRY_RUN" -ne 1 ]]; then
  mkdir -p "$ABS_TRANSFER_DIR"
fi

MANIFEST="${ABS_TRANSFER_DIR}/staging_manifest.tsv"
TMP_MANIFEST=""
cleanup() {
  if [[ -n "$TMP_MANIFEST" && -f "$TMP_MANIFEST" ]]; then
    rm -f "$TMP_MANIFEST"
  fi
}
trap cleanup EXIT

if [[ "$DRY_RUN" -ne 1 ]]; then
  TMP_MANIFEST="$(mktemp "${ABS_TRANSFER_DIR}/staging_manifest.XXXXXX")"
  printf 'feature\tgroup\tsource\tstaged_file\n' > "$TMP_MANIFEST"
fi

for idx in "${!PLAN_SRCS[@]}"; do
  src="${PLAN_SRCS[$idx]}"
  dest="${PLAN_DESTS[$idx]}"
  feature="${PLAN_FEATURES[$idx]}"
  group="${PLAN_GROUPS[$idx]}"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf 'DRY-RUN\t%s\t%s\n' "$src" "$dest"
  else
    mkdir -p "$(dirname -- "$dest")"
    echo "staging $src"
    rsync -aL -- "$src" "$dest"
    printf '%s\t%s\t%s\t%s\n' "$feature" "$group" "$src" "$dest" >> "$TMP_MANIFEST"
  fi
done

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry run complete. Files identified: ${#PLAN_SRCS[@]}"
else
  mv "$TMP_MANIFEST" "$MANIFEST"
  TMP_MANIFEST=""
  echo "Staged files: ${#PLAN_SRCS[@]}"
  echo "Transfer directory: $ABS_TRANSFER_DIR"
  echo "Manifest: $MANIFEST"
fi

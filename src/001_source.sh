#!/bin/bash
# Script: src/001_source.sh
# Purpose: Set up the lightweight Nextflow environment and bundled pipeline tools.

set -euo pipefail

SCRIPT_SOURCE="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(dirname -- "$SCRIPT_SOURCE")" && pwd)"
DEFAULT_PROJ_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJ_ROOT="${GWAS_PROJECT_ROOT:-$DEFAULT_PROJ_ROOT}"
ENV_NAME="nf_master"

source_env_file() {
  local env_file="$1"
  if [[ -f "$env_file" ]]; then
    set -a
    source "$env_file"
    set +a
  fi
}

source_env_file "${PWD}/.env"
if [[ "$DEFAULT_PROJ_ROOT" != "$PWD" ]]; then
  source_env_file "${DEFAULT_PROJ_ROOT}/.env"
fi

PIPELINE_DIR="${PROJ_ROOT}/pipeline"
TOOLS_DIR="${PIPELINE_DIR}/tools"
METAL_URL="https://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz"
FORCE=0
SETUP_ENV=1
SETUP_TOOLS=1

choose_conda_cmd() {
  if command -v mamba >/dev/null 2>&1; then
    printf '%s\n' "mamba"
    return 0
  fi

  if command -v conda >/dev/null 2>&1; then
    printf '%s\n' "conda"
    return 0
  fi

  echo "ERROR: Neither conda nor mamba is available in your PATH." >&2
  exit 1
}

setup_nextflow_env() {
  local conda_cmd
  conda_cmd="$(choose_conda_cmd)"

  echo "=========================================="
  echo " Setting up Nextflow Master Environment"
  echo "=========================================="

  if "$conda_cmd" env list | awk 'NR > 2 {print $1}' | grep -Fxq "$ENV_NAME"; then
    echo "Environment '${ENV_NAME}' already exists. Skipping creation."
    return 0
  fi

  echo "Creating environment '${ENV_NAME}' with ${conda_cmd}..."
  "$conda_cmd" create -n "$ENV_NAME" -c conda-forge -c bioconda nextflow openjdk=21 python=3.9 -y

  echo "=========================================="
  echo "Environment ready. Activate with: 'conda activate ${ENV_NAME}'"
  echo "=========================================="
}

install_metal() {
  echo "=========================================="
  echo " Installing Bundled METAL Binary"
  echo "=========================================="

  mkdir -p "$TOOLS_DIR"
  ARCHIVE="${TOOLS_DIR}/Linux-metal.tar.gz"

  if [[ -f "$ARCHIVE" && "$FORCE" -ne 1 ]]; then
    echo "Archive already exists: $ARCHIVE (use --force to re-download)"
  else
    echo "Downloading METAL from: $METAL_URL"
    if command -v curl >/dev/null 2>&1; then
      curl -L -o "$ARCHIVE" "$METAL_URL"
    elif command -v wget >/dev/null 2>&1; then
      wget -O "$ARCHIVE" "$METAL_URL"
    else
      echo "ERROR: curl or wget is required to download METAL."
      exit 1
    fi
  fi

  TMPDIR="$(mktemp -d -p "$TOOLS_DIR")"
  tar -xzf "$ARCHIVE" -C "$TMPDIR"

  METAL_BIN="$(find "$TMPDIR" -type f -name metal -perm -u+x | head -n 1)"
  if [[ -z "$METAL_BIN" ]]; then
    METAL_BIN="$(find "$TMPDIR" -type f -name metal | head -n 1)"
  fi

  if [[ -z "$METAL_BIN" ]]; then
    echo "ERROR: Could not locate METAL binary in extracted archive."
    rm -rf "$TMPDIR"
    exit 1
  fi

  mkdir -p "${TOOLS_DIR}/metal"
  cp "$METAL_BIN" "${TOOLS_DIR}/metal/metal"
  chmod +x "${TOOLS_DIR}/metal/metal"

  rm -rf "$TMPDIR"

  echo "METAL installed to: ${TOOLS_DIR}/metal/metal"
}

while [[ $# -gt 0 ]]; do
  case $1 in
    --env-name) ENV_NAME="$2"; shift 2 ;;
    --tools-dir) TOOLS_DIR="$2"; shift 2 ;;
    --url) METAL_URL="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --env-only)
      SETUP_TOOLS=0
      shift
      ;;
    --tools-only)
      SETUP_ENV=0
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [--env-name NAME] [--tools-dir PATH] [--url URL] [--force] [--env-only|--tools-only]"
      echo "Default env name: ${ENV_NAME}"
      echo "Default install location: ${TOOLS_DIR}"
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ "$SETUP_ENV" -eq 1 ]]; then
  setup_nextflow_env
fi

if [[ "$SETUP_TOOLS" -eq 1 ]]; then
  install_metal
fi

#!/bin/bash
# Download and install METAL binary into tools/

set -euo pipefail

# Robust environment sourcing
# Try .env in current directory, then relative to script, then using GWAS_PROJECT_ROOT
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

# Set project root to current directory if not set via .env
PROJ_ROOT="${GWAS_PROJECT_ROOT:-$(pwd)}"
TOOLS_DIR="${PROJ_ROOT}/tools"
METAL_URL="https://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz"
FORCE=0

while [[ $# -gt 0 ]]; do
  case $1 in
    --tools-dir) TOOLS_DIR="$2"; shift 2 ;;
    --url) METAL_URL="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help)
      echo "Usage: $0 [--tools-dir PATH] [--url URL] [--force]"
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

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

#!/bin/bash
# Standalone METAL test for cases vs controls within a study/protein.
# This script does NOT call any other pipeline scripts.

set -euo pipefail

OUTDIR="/data/Epic/subprojects/Somalogic/work/GWAS/analysis/test3"
STUDY="Neuro_01"
PROTEIN="10620-21"
METAL_BIN="metal"
CASES_FILE=""
CONTROLS_FILE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --study) STUDY="$2"; shift 2 ;;
    --protein) PROTEIN="$2"; shift 2 ;;
    --metal) METAL_BIN="$2"; shift 2 ;;
    --cases) CASES_FILE="$2"; shift 2 ;;
    --controls) CONTROLS_FILE="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--outdir PATH] [--study ID] [--protein ID] [--metal PATH] [--cases FILE] [--controls FILE]"
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ -z "$CASES_FILE" ]]; then
  CASES_FILE="${OUTDIR}/${STUDY}/${PROTEIN}/cases/009_QC/gwas.tsv.gz"
fi
if [[ -z "$CONTROLS_FILE" ]]; then
  CONTROLS_FILE="${OUTDIR}/${STUDY}/${PROTEIN}/controls/009_QC/gwas.tsv.gz"
fi

META_DIR="${OUTDIR}/${STUDY}/${PROTEIN}/meta"

if [[ ! -f "$CASES_FILE" || ! -f "$CONTROLS_FILE" ]]; then
  echo "ERROR: Missing case/control GWAS files."
  echo "Expected:"
  echo "  $CASES_FILE"
  echo "  $CONTROLS_FILE"
  exit 1
fi

mkdir -p "$META_DIR"
TMPDIR=$(mktemp -d -p "$META_DIR")
INPUT_LIST="${TMPDIR}/gwama_input.txt"

format_for_metal () {
  local infile="$1"
  local outfile="$2"

  zcat "$infile" | awk -v FS='[ \t]+' -v OFS='\t' '
    NR==1 {
      sub(/^\xef\xbb\xbf/, "", $1)
      for (i=1; i<=NF; i++) {
        gsub(/\r/, "", $i)
        key = toupper($i)
        h[key]=i
      }
      id = (h["ID"] ? h["ID"] : (h["MARKERNAME"] ? h["MARKERNAME"] : 0))
      ea = (h["EA"] ? h["EA"] : (h["ALLELE1"] ? h["ALLELE1"] : 0))
      nea = (h["OA"] ? h["OA"] : (h["ALLELE0"] ? h["ALLELE0"] : 0))
      eaf = (h["EAF"] ? h["EAF"] : (h["A1FREQ"] ? h["A1FREQ"] : 0))
      beta = (h["BETA"] ? h["BETA"] : 0)
      se = (h["SE"] ? h["SE"] : 0)
      n = (h["N"] ? h["N"] : 0)

      if (id==0 && $1=="CHR" && $2=="POS" && $3=="ID") {
        id=3; ea=4; nea=5; eaf=6; beta=7; se=8; n=10
      } else if (id==0 && $1=="CHROM" && $2=="GENPOS" && $3=="ID") {
        id=3; ea=5; nea=4; eaf=6; beta=10; se=11; n=8
      }

      if (id==0 || ea==0 || nea==0 || eaf==0 || beta==0 || se==0 || n==0) {
        print "ERROR: Required columns missing for METAL in " FILENAME > "/dev/stderr"
        exit 1
      }
      print "MARKERNAME\tEA\tNEA\tEAF\tBETA\tSE\tN"
      next
    }
    NR>1 {
      print $(id), $(ea), $(nea), $(eaf), $(beta), $(se), $(n)
    }
  ' > "$outfile"
} 

format_for_metal "$CASES_FILE" "${TMPDIR}/input_0_gwas.tsv.txt"
format_for_metal "$CONTROLS_FILE" "${TMPDIR}/input_1_gwas.tsv.txt"

printf "%s\n" "${TMPDIR}/input_0_gwas.tsv.txt" "${TMPDIR}/input_1_gwas.tsv.txt" > "$INPUT_LIST"

metal_script="${TMPDIR}/metal_script.txt"
cat > "$metal_script" <<EOF
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
MARKER ID
ALLELE EA OA
EFFECT BETA
STDERR SE
PVALUE P
WEIGHT N
FREQ EAF
EOF

while read -r f; do
  echo "PROCESS $f" >> "$metal_script"
done < "$INPUT_LIST"

cat >> "$metal_script" <<EOF
OUTFILE ${META_DIR}/meta .tbl
ANALYZE
QUIT
EOF

"$METAL_BIN" "$metal_script" > "${META_DIR}/meta.log" 2>&1

if [[ -f "${META_DIR}/meta.tbl" ]]; then
  gzip -f "${META_DIR}/meta.tbl"
elif [[ -f "${META_DIR}/meta1.tbl" ]]; then
  gzip -f "${META_DIR}/meta1.tbl"
fi

rm -rf "$TMPDIR"

echo "METAL meta completed:"
echo "  ${META_DIR}/meta.tbl.gz (or meta1.tbl.gz)"
echo "  ${META_DIR}/meta.log"

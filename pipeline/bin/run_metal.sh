#!/bin/bash
set -euo pipefail

input_files=""
protein_id=""
study=""
group=""
metal_binary="metal"
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --input_files) input_files="$2"; shift 2 ;;
    --protein_id) protein_id="$2"; shift 2 ;;
    --study) study="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --metal_binary) metal_binary="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"
tmpdir=$(mktemp -d -p "$outdir")
cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
metal_qc_script="${script_dir}/metal_qc.py"

if [[ ! -x "$metal_binary" ]]; then
  echo "ERROR: METAL binary not found or not executable: $metal_binary" > "${outdir}/meta.log"
  exit 1
fi

if [[ ! -f "$metal_qc_script" ]]; then
  echo "ERROR: METAL QC helper not found: $metal_qc_script" > "${outdir}/meta.log"
  exit 1
fi

IFS=',' read -r -a FILES <<< "$input_files"
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "ERROR: No input files were provided to METAL" > "${outdir}/meta.log"
  exit 1
fi

input_list="${tmpdir}/metal_input.txt"
standardized_header=$'CHR\tPOS\tID\tEA\tOA\tEAF\tBETA\tSE\tP\tN'

read_header() {
  local infile="$1"

  if gzip -t "$infile" >/dev/null 2>&1; then
    python3 - "$infile" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "rt", encoding="utf-8", errors="replace") as handle:
    print(handle.readline().rstrip("\r\n"))
PY
  else
    sed -n '1p' "$infile" | tr -d '\r'
  fi
}

format_standardized_input() {
  local infile="$1"
  local outfile="$2"
  local compressed="$3"

  if [[ "$compressed" == "1" ]]; then
    gzip -cd -- "$infile" | awk 'BEGIN{FS=OFS="\t"}
      NR==1 { print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"; next }
      NF == 0 { next }
      {
        chr = $1
        gsub(/^[Cc][Hh][Rr]/, "", chr)
        print $3,$4,$5,$6,$7,$8,$9,$10,chr,$2
      }' > "$outfile"
  else
    awk 'BEGIN{FS=OFS="\t"}
      NR==1 { print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"; next }
      NF == 0 { next }
      {
        chr = $1
        gsub(/^[Cc][Hh][Rr]/, "", chr)
        print $3,$4,$5,$6,$7,$8,$9,$10,chr,$2
      }' "$infile" > "$outfile"
  fi
}

format_generic_input() {
  local infile="$1"
  local outfile="$2"
  local compressed="$3"

  if [[ "$compressed" == "1" ]]; then
    gzip -cd -- "$infile" | awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
      {
        if ($0 ~ /^#/ || NF == 0) next
        if (!header_done) {
          for (i = 1; i <= NF; i++) {
            gsub(/"/, "", $i); gsub(/\r/, "", $i)
            t = toupper($i); gsub(/[^A-Z0-9]/, "", t)
            if (t == "ID" || t == "MARKERNAME" || t == "SNP") id = i
            else if (t == "EA" || t == "ALLELE1" || t == "A1" || t == "EFFECTALLELE") ea = i
            else if (t == "OA" || t == "ALLELE2" || t == "A2" || t == "OTHERALLELE") oa = i
            else if (t == "EAF" || t == "FREQ1" || t == "A1FREQ" || t == "FREQ" || t == "FREQA1") eaf = i
            else if (t == "BETA" || t == "EFFECT") beta = i
            else if (t == "SE" || t == "STDERR" || t == "STANDARDERROR") se = i
            else if (t == "P" || t == "PVALUE" || t == "PVAL") p = i
            else if (t == "N" || t == "WEIGHT" || t == "TOTALSAMPLESIZE" || t == "NEFF") n = i
            else if (t == "CHR" || t == "CHROM") chr = i
            else if (t == "POS" || t == "POSITION" || t == "GENPOS") pos = i
          }
          if (id == 0 || ea == 0 || oa == 0 || eaf == 0 || beta == 0 || se == 0 || p == 0) {
            print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
            exit 2
          }
          print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"
          header_done = 1
          next
        }
        nval = (n > 0 ? $(n) : "1")
        chrval = (chr > 0 ? $(chr) : "")
        posval = (pos > 0 ? $(pos) : "")
        gsub(/^[Cc][Hh][Rr]/, "", chrval)
        print $(id),$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval,chrval,posval
      }' > "$outfile"
  else
    awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
      {
        if ($0 ~ /^#/ || NF == 0) next
        if (!header_done) {
          for (i = 1; i <= NF; i++) {
            gsub(/"/, "", $i); gsub(/\r/, "", $i)
            t = toupper($i); gsub(/[^A-Z0-9]/, "", t)
            if (t == "ID" || t == "MARKERNAME" || t == "SNP") id = i
            else if (t == "EA" || t == "ALLELE1" || t == "A1" || t == "EFFECTALLELE") ea = i
            else if (t == "OA" || t == "ALLELE2" || t == "A2" || t == "OTHERALLELE") oa = i
            else if (t == "EAF" || t == "FREQ1" || t == "A1FREQ" || t == "FREQ" || t == "FREQA1") eaf = i
            else if (t == "BETA" || t == "EFFECT") beta = i
            else if (t == "SE" || t == "STDERR" || t == "STANDARDERROR") se = i
            else if (t == "P" || t == "PVALUE" || t == "PVAL") p = i
            else if (t == "N" || t == "WEIGHT" || t == "TOTALSAMPLESIZE" || t == "NEFF") n = i
            else if (t == "CHR" || t == "CHROM") chr = i
            else if (t == "POS" || t == "POSITION" || t == "GENPOS") pos = i
          }
          if (id == 0 || ea == 0 || oa == 0 || eaf == 0 || beta == 0 || se == 0 || p == 0) {
            print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
            exit 2
          }
          print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"
          header_done = 1
          next
        }
        nval = (n > 0 ? $(n) : "1")
        chrval = (chr > 0 ? $(chr) : "")
        posval = (pos > 0 ? $(pos) : "")
        gsub(/^[Cc][Hh][Rr]/, "", chrval)
        print $(id),$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval,chrval,posval
      }' "$infile" > "$outfile"
  fi
}

for i in "${!FILES[@]}"; do
  f="${FILES[$i]}"
  bn=$(basename "$f")
  out_txt="${tmpdir}/input_${i}_${bn%.gz}.txt"
  compressed="0"
  header=$(read_header "$f")

  if gzip -t "$f" >/dev/null 2>&1; then
    compressed="1"
  fi

  if [[ "$header" == "${standardized_header}"* ]]; then
    format_standardized_input "$f" "$out_txt" "$compressed"
  else
    format_generic_input "$f" "$out_txt" "$compressed"
  fi

  echo "$out_txt" >> "$input_list"
done

metal_script="${tmpdir}/metal_script.txt"
cat > "$metal_script" <<EOF
# TRACKPOSITIONS ON ensures METAL correctly tracks and outputs Chromosome/Position
TRACKPOSITIONS ON
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N
MARKER ID
ALLELE EA OA
EFFECT BETA
STDERR SE
PVALUE P
WEIGHT N
FREQ EAF
CHROMOSOMELABEL CHR
POSITIONLABEL POS
EOF

while read -r f; do
  echo "PROCESS $f" >> "$metal_script"
done < "$input_list"

cat >> "$metal_script" <<EOF
OUTFILE ${outdir}/meta .tbl
ANALYZE
QUIT
EOF

"$metal_binary" "$metal_script" > "${outdir}/meta.log" 2>&1

tbl_file=""
if [[ -f "${outdir}/meta.tbl" ]]; then
  tbl_file="${outdir}/meta.tbl"
else
  tbl_file=$(ls -1 "${outdir}"/meta*.tbl 2>/dev/null | head -n 1 || true)
fi

if [[ -z "$tbl_file" ]]; then
  echo "ERROR: METAL did not produce a .tbl output for ${protein_id}" >> "${outdir}/meta.log"
  exit 1
fi

stem="${tbl_file%.tbl}"
if [[ -f "${stem}.tbl.info" && "${stem}.tbl.info" != "${outdir}/meta.tbl.info" ]]; then
  mv -f "${stem}.tbl.info" "${outdir}/meta.tbl.info"
fi

python3 "$metal_qc_script" \
  --input_file "$tbl_file" \
  --protein_id "$protein_id" \
  --study "$study" \
  --group "$group" \
  --outdir "$outdir"

rm -f "$tbl_file"

echo "METAL for $protein_id ($group) complete."

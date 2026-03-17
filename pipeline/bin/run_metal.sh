#!/bin/bash
set -euo pipefail

input_files=""
protein_id=""
group=""
metal_binary="metal"
outdir=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --input_files) input_files="$2"; shift 2 ;;
    --protein_id) protein_id="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --metal_binary) metal_binary="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"
tmpdir=$(mktemp -d -p "$outdir")

if [[ ! -x "$metal_binary" ]]; then
  echo "ERROR: METAL binary not found or not executable: $metal_binary" > "${outdir}/meta.log"
  exit 1
fi

IFS=',' read -ra FILES <<< "$input_files"
input_list="${tmpdir}/metal_input.txt"

format_for_metal () {
  local infile="$1"
  local outfile="$2"

  if [[ "$infile" == *.gz ]]; then
    zcat "$infile" | awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
      {
        if ($0 ~ /^#/ || NF==0) next
        if (!header_done) {
          for(i=1;i<=NF;i++){
            gsub(/"/,"",$i)
            gsub(/\r/,"",$i)
            t=toupper($i); gsub(/[^A-Z0-9]/,"",t)
            if(t=="ID"||t=="MARKERNAME"||t=="SNP") id=i
            else if(t=="EA"||t=="ALLELE1"||t=="A1"||t=="EFFECTALLELE") ea=i
            else if(t=="OA"||t=="ALLELE2"||t=="A2"||t=="OTHERALLELE") oa=i
            else if(t=="EAF"||t=="FREQ1"||t=="A1FREQ"||t=="FREQ"||t=="FREQA1") eaf=i
            else if(t=="BETA"||t=="EFFECT") beta=i
            else if(t=="SE"||t=="STDERR"||t=="STANDARDERROR") se=i
            else if(t=="P"||t=="PVALUE"||t=="PVAL") p=i
            else if(t=="N"||t=="WEIGHT"||t=="TOTALSAMPLESIZE"||t=="NEFF") n=i
          }
          if(id==0||ea==0||oa==0||eaf==0||beta==0||se==0||p==0){
            print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
            exit 2
          }
          print "ID","EA","OA","EAF","BETA","SE","P","N"
          header_done=1
          next
        }
        nval = (n>0 ? $(n) : "1")
        print $(id),$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval
      }
      END{
        if(!header_done){
          print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
          exit 2
        }
      }' > "$outfile"
  else
    awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
      {
        if ($0 ~ /^#/ || NF==0) next
        if (!header_done) {
          for(i=1;i<=NF;i++){
            gsub(/"/,"",$i)
            gsub(/\r/,"",$i)
            t=toupper($i); gsub(/[^A-Z0-9]/,"",t)
            if(t=="ID"||t=="MARKERNAME"||t=="SNP") id=i
            else if(t=="EA"||t=="ALLELE1"||t=="A1"||t=="EFFECTALLELE") ea=i
            else if(t=="OA"||t=="ALLELE2"||t=="A2"||t=="OTHERALLELE") oa=i
            else if(t=="EAF"||t=="FREQ1"||t=="A1FREQ"||t=="FREQ"||t=="FREQA1") eaf=i
            else if(t=="BETA"||t=="EFFECT") beta=i
            else if(t=="SE"||t=="STDERR"||t=="STANDARDERROR") se=i
            else if(t=="P"||t=="PVALUE"||t=="PVAL") p=i
            else if(t=="N"||t=="WEIGHT"||t=="TOTALSAMPLESIZE"||t=="NEFF") n=i
          }
          if(id==0||ea==0||oa==0||eaf==0||beta==0||se==0||p==0){
            print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
            exit 2
          }
          print "ID","EA","OA","EAF","BETA","SE","P","N"
          header_done=1
          next
        }
        nval = (n>0 ? $(n) : "1")
        print $(id),$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval
      }
      END{
        if(!header_done){
          print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
          exit 2
        }
      }' "$infile" > "$outfile"
  fi
}

# Prepare plain-text files for METAL
for i in "${!FILES[@]}"; do
  f="${FILES[$i]}"
  bn=$(basename "$f")
  out_txt="${tmpdir}/input_${i}_${bn%.gz}.txt"
  format_for_metal "$f" "$out_txt"
  echo "$out_txt" >> "$input_list"
done

metal_script="${tmpdir}/metal_script.txt"
cat > "$metal_script" <<EOF
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
  tbl_file=$(ls -1 "${outdir}"/meta*.tbl 2>/dev/null | head -n1 || true)
fi

if [[ -n "$tbl_file" ]]; then
  stem="${tbl_file%.tbl}"
  info_file="${stem}.tbl.info"

  if [[ "$tbl_file" != "${outdir}/meta.tbl" ]]; then
    mv -f "$tbl_file" "${outdir}/meta.tbl"
  fi
  if [[ -f "$info_file" ]]; then
    mv -f "$info_file" "${outdir}/meta.tbl.info"
  fi

  # Standardize output columns so METAL outputs can be re-used as inputs
  # Log first non-empty line to aid debugging
  awk 'NF>0 && $0 !~ /^#/ {print "META_HEADER_LINE: " $0; exit}' "${outdir}/meta.tbl" >> "${outdir}/meta.log" || true

  if awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
  {
    if ($0 ~ /^#/ || NF==0) next
    if (!header_done) {
      for(i=1;i<=NF;i++){
        gsub(/"/,"",$i)
        gsub(/\r/,"",$i)
        t=toupper($i); gsub(/[^A-Z0-9]/,"",t)
        if(t=="MARKERNAME"||t=="ID"||t=="SNP") id=i
        else if(t=="ALLELE1"||t=="EA"||t=="A1"||t=="EFFECTALLELE") ea=i
        else if(t=="ALLELE2"||t=="OA"||t=="A2"||t=="OTHERALLELE") oa=i
        else if(t=="FREQ1"||t=="EAF"||t=="A1FREQ"||t=="FREQ"||t=="FREQA1") eaf=i
        else if(t=="EFFECT"||t=="BETA") beta=i
        else if(t=="STDERR"||t=="SE"||t=="STANDARDERROR") se=i
        else if(t=="PVALUE"||t=="P"||t=="PVAL") p=i
        else if(t=="WEIGHT"||t=="N"||t=="TOTALSAMPLESIZE"||t=="NEFF") n=i
      }
      if(id==0||ea==0||oa==0||eaf==0||beta==0||se==0||p==0){
        print "ERROR: Required columns missing in METAL output" > "/dev/stderr"
        exit 2
      }
      print "ID","EA","OA","EAF","BETA","SE","P","N"
      header_done=1
      next
    }
    nval = (n>0 ? $(n) : "1")
    print $(id),$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval
  }
  END{
    if(!header_done){
      print "ERROR: Required columns missing in METAL output" > "/dev/stderr"
      exit 2
    }
  }' \
  "${outdir}/meta.tbl" > "${outdir}/meta.tbl.std"; then
    mv -f "${outdir}/meta.tbl.std" "${outdir}/meta.tbl"
  else
    echo "WARN: Could not standardize METAL output columns; keeping original meta.tbl" >> "${outdir}/meta.log"
    rm -f "${outdir}/meta.tbl.std"
  fi

  if [[ ! -f "${outdir}/meta.tbl.info" ]]; then
    echo "WARN: meta.tbl.info missing" >> "${outdir}/meta.log"
    touch "${outdir}/meta.tbl.info"
  fi

  gzip -f "${outdir}/meta.tbl"
fi

if [[ ! -f "${outdir}/meta.tbl.gz" ]]; then
  echo "ERROR: METAL did not produce expected output files. Log:" >> "${outdir}/meta.log"
  cat "${outdir}/meta.log"
  exit 1
fi

# Cleanup: keep only meta.tbl.gz, meta.tbl.info, meta.log
for f in "${outdir}"/*; do
  base="$(basename "$f")"
  case "$base" in
    meta.tbl.gz|meta.tbl.info|meta.log) ;;
    *)
      if [[ -d "$f" ]]; then
        rm -rf "$f"
      else
        rm -f "$f"
      fi
      ;;
  esac
done

rm -rf "$tmpdir"

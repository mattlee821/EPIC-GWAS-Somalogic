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

# Pre-process inputs to ensure clean headers for METAL
format_for_metal () {
  local infile="$1"
  local outfile="$2"

  if gzip -t "$infile" >/dev/null 2>&1; then
    gzip -cd "$infile" | awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
      {
        if ($0 ~ /^#/ || NF==0) next
        if (!header_done) {
          for(i=1;i<=NF;i++){
            gsub(/"/,"",$i); gsub(/\r/,"",$i)
            t=toupper($i); gsub(/[^A-Z0-9]/,"",t)
            if(t=="ID"||t=="MARKERNAME"||t=="SNP") id=i
            else if(t=="EA"||t=="ALLELE1"||t=="A1"||t=="EFFECTALLELE") ea=i
            else if(t=="OA"||t=="ALLELE2"||t=="A2"||t=="OTHERALLELE") oa=i
            else if(t=="EAF"||t=="FREQ1"||t=="A1FREQ"||t=="FREQ"||t=="FREQA1") eaf=i
            else if(t=="BETA"||t=="EFFECT") beta=i
            else if(t=="SE"||t=="STDERR"||t=="STANDARDERROR") se=i
            else if(t=="P"||t=="PVALUE"||t=="PVAL") p=i
            else if(t=="N"||t=="WEIGHT"||t=="TOTALSAMPLESIZE"||t=="NEFF") n=i
            else if(t=="CHR"||t=="CHROM") chr=i
            else if(t=="POS"||t=="POSITION"||t=="GENPOS") pos=i
          }
          if(id==0||ea==0||oa==0||eaf==0||beta==0||se==0||p==0){
            print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
            exit 2
          }
          print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"
          header_done=1
          next
        }
        nval = (n>0 ? $(n) : "1")
        chrval = (chr>0 ? $(chr) : "")
        posval = (pos>0 ? $(pos) : "")
        # Remove "chr" prefix if present
        gsub(/^[Cc][Hh][Rr]/, "", chrval)
        print $(id),$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval,chrval,posval
      }' > "$outfile"
  else
    awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
      {
        if ($0 ~ /^#/ || NF==0) next
        if (!header_done) {
          for(i=1;i<=NF;i++){
            gsub(/"/,"",$i); gsub(/\r/,"",$i)
            t=toupper($i); gsub(/[^A-Z0-9]/,"",t)
            if(t=="ID"||t=="MARKERNAME"||t=="SNP") id=i
            else if(t=="EA"||t=="ALLELE1"||t=="A1"||t=="EFFECTALLELE") ea=i
            else if(t=="OA"||t=="ALLELE2"||t=="A2"||t=="OTHERALLELE") oa=i
            else if(t=="EAF"||t=="FREQ1"||t=="A1FREQ"||t=="FREQ"||t=="FREQA1") eaf=i
            else if(t=="BETA"||t=="EFFECT") beta=i
            else if(t=="SE"||t=="STDERR"||t=="STANDARDERROR") se=i
            else if(t=="P"||t=="PVALUE"||t=="PVAL") p=i
            else if(t=="N"||t=="WEIGHT"||t=="TOTALSAMPLESIZE"||t=="NEFF") n=i
            else if(t=="CHR"||t=="CHROM") chr=i
            else if(t=="POS"||t=="POSITION"||t=="GENPOS") pos=i
          }
          if(id==0||ea==0||oa==0||eaf==0||beta==0||se==0||p==0){
            print "ERROR: Required columns missing in input for METAL" > "/dev/stderr"
            exit 2
          }
          print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"
          header_done=1
          next
        }
        nval = (n>0 ? $(n) : "1")
        chrval = (chr>0 ? $(chr) : "")
        posval = (pos>0 ? $(pos) : "")
        gsub(/^[Cc][Hh][Rr]/, "", chrval)
        print $(id),$(ea),$(oa),$(eaf),$(beta),$(se),$(p),nval,chrval,posval
      }' "$infile" > "$outfile"
  fi
}

for i in "${!FILES[@]}"; do
  f="${FILES[$i]}"
  bn=$(basename "$f")
  out_txt="${tmpdir}/input_${i}_${bn%.gz}.txt"
  format_for_metal "$f" "$out_txt"
  echo "$out_txt" >> "$input_list"
done

metal_script="${tmpdir}/metal_script.txt"
cat > "$metal_script" <<EOF
# TRACKPOSITIONS ON ensures METAL correctly tracks and outputs Chromosome/Position
TRACKPOSITIONS ON
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
# TotalSampleSize setup
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N
# Required mappings
MARKER ID
ALLELE EA OA
EFFECT BETA
STDERR SE
PVALUE P
WEIGHT N
FREQ EAF
# Native Coordinate Tracking
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
  tbl_file=$(ls -1 "${outdir}"/meta*.tbl 2>/dev/null | head -n1 || true)
fi

if [[ -n "$tbl_file" ]]; then
  stem="${tbl_file%.tbl}"
  # Ensure the .info file is correctly named
  if [[ -f "${stem}.tbl.info" && "${stem}.tbl.info" != "${outdir}/meta.tbl.info" ]]; then
      mv -f "${stem}.tbl.info" "${outdir}/meta.tbl.info"
  fi

  # Final Standardization Step (Mapping Chromosome/Position columns)
  awk 'BEGIN{OFS="\t"; FS="[ \t]+"}
  {
    if ($0 ~ /^#/ || NF==0) next
    if (!header_done) {
      for(i=1;i<=NF;i++){
        gsub(/"/,"",$i); gsub(/\r/,"",$i)
        t=toupper($i); gsub(/[^A-Z0-9]/,"",t)
        if(t=="MARKERNAME"||t=="ID"||t=="SNP") id_idx=i
        else if(t=="ALLELE1"||t=="EA"||t=="A1"||t=="EFFECTALLELE") ea_idx=i
        else if(t=="ALLELE2"||t=="OA"||t=="A2"||t=="OTHERALLELE") oa_idx=i
        else if(t=="FREQ1"||t=="EAF"||t=="A1FREQ"||t=="FREQ"||t=="FREQA1") eaf_idx=i
        else if(t=="EFFECT"||t=="BETA") beta_idx=i
        else if(t=="STDERR"||t=="SE"||t=="STANDARDERROR") se_idx=i
        else if(t=="PVALUE"||t=="P"||t=="PVAL") p_idx=i
        else if(t=="WEIGHT"||t=="N"||t=="TOTALSAMPLESIZE"||t=="NEFF") n_idx=i
        else if(t=="CHROMOSOME"||t=="CHR") chr_idx=i
        else if(t=="POSITION"||t=="POS") pos_idx=i
      }
      print "ID","EA","OA","EAF","BETA","SE","P","N","CHR","POS"
      header_done=1
      next
    }
    # Fetch result variables
    nval = (n_idx>0 ? $(n_idx) : "1")
    chrval = (chr_idx>0 ? $(chr_idx) : "")
    posval = (pos_idx>0 ? $(pos_idx) : "")
    
    # Fallback to string split if somehow coordinates are missing
    if(chrval=="" || posval==""){
        split($(id_idx), a, /[_:]/)
        if(chrval=="") chrval=a[1]
        if(posval=="") posval=a[2]
    }
    gsub(/^[Cc][Hh][Rr]/, "", chrval)
    print $(id_idx),$(ea_idx),$(oa_idx),$(eaf_idx),$(beta_idx),$(se_idx),$(p_idx),nval,chrval,posval
  }' "$tbl_file" > "${outdir}/meta.tbl.std"
  
  mv -f "${outdir}/meta.tbl.std" "${outdir}/meta.tbl"
  gzip -f "${outdir}/meta.tbl"
fi

# Cleanup
rm -rf "$tmpdir"
echo "METAL for $protein_id ($group) complete."

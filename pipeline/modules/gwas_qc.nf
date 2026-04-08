process GWAS_QC {
  label        'low'
  conda        "${projectDir}/envs/py_gwas.yml"
  publishDir   params.outdir, mode: 'copy', saveAs: { name ->
    def parts = name.tokenize('/')
    def protein_id = parts ? parts[0] : "UNKNOWN"
    def base = parts ? parts[-1] : name
    "${study_id}/${protein_id}/${group}/009_QC/${base}"
  }

  input:
  tuple val(study_id), val(group), val(protein_ids), path(regenie_files), val(mapping)

  output:
  tuple val(study_id), val(group), val(protein_ids), path("*/gwas.tsv.gz"), path("*/metrics.tsv"), emit: all_out
  path "*/gwas.tsv.gz", emit: gwas_files

  path "*/metrics.tsv", emit: metrics_files

  script:
  """
  echo -e "${mapping}" > mapping.txt
  while IFS=':' read -r prot files; do
    mkdir -p "\${prot}"
    python ${projectDir}/bin/gwas_qc.py \\
      --input_files "\${files}" \\
      --protein_id "\${prot}" \\
      --study "${study_id}" \\
      --group "${group}" \\
      --info_score ${params.info_score} \\
      --outdir "\${prot}"
  done < mapping.txt
  """



}

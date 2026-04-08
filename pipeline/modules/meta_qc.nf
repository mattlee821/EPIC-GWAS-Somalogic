process META_QC {
  label        'medium'
  conda        "${projectDir}/envs/py_gwas.yml"
  publishDir   {
    def pid = (protein_ids instanceof java.util.Collection && protein_ids) ? protein_ids[0] : protein_ids
    if (study_id == "meta-study" || study_id == "meta_study") {
      "${params.outdir}/meta/${pid}/${analysis}/013_QC"
    } else {
      "${params.outdir}/${study_id}/${pid}/${analysis}/011_QC"
    }
  }, mode: 'copy'

  input:
  tuple val(study_id), val(analysis), val(protein_ids), path(meta_tbls), val(mapping)

  output:
  tuple val(study_id), val(analysis), val(protein_ids), path("meta.tsv.gz"), path("metrics.tsv"), emit: all_out
  path "meta.tsv.gz", emit: gwas_files
  path "metrics.tsv", emit: metrics_files

  script:
  """
  echo -e "${mapping}" > mapping.txt
  while IFS=':' read -r prot file; do
    python ${projectDir}/bin/metal_qc.py \\
      --input_file "\${file}" \\
      --protein_id "\${prot}" \\
      --study "${study_id}" \\
      --group "${analysis}" \\
      --outdir .
  done < mapping.txt
  """


}

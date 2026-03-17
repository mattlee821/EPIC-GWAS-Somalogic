process META_QC {
  label        'medium'
  conda        "${projectDir}/envs/py_gwas.yml"
  publishDir   { study_id == "meta-study" ? "${params.outdir}/meta/${protein_id}/${analysis}/QC" : "${params.outdir}/${study_id}/${protein_id}/meta/QC" }, mode: 'copy'

  input:
  tuple val(study_id), val(analysis), val(protein_id), path(meta_tbl)

  output:
  tuple val(study_id), val(analysis), val(protein_id), path("gwas.tsv.gz"), path("metrics.tsv"), path("figure.png")

  script:
  """
  python ${projectDir}/bin/finalise.py \\
    --input_files "${meta_tbl}" \\
    --study ${study_id} \\
    --group ${analysis} \\
    --protein_id ${protein_id} \\
    --info_score ${params.info_score} \\
    --outdir .
  """
}

process POST_GWAS_QC {
  label        'low'
  conda        "${projectDir}/envs/py_gwas.yml"
  publishDir   { "${params.outdir}/${study_id}/${protein_id}/${group}/009_QC" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), val(protein_id), path(regenie_files)

  output:
  tuple val(study_id), val(group), val(protein_id),
        path("gwas.tsv.gz"),
        path("metrics.tsv"),
        path("figure.png")

  script:
  def files_arg = regenie_files.join(',')
  """
  # 1. Run the python script
  python ${projectDir}/bin/finalise.py \\
    --input_files "${files_arg}" \\
    --study ${study_id} \\
    --group ${group} \\
    --protein_id ${protein_id} \\
    --info_score ${params.info_score} \\
    --outdir .

  """
}

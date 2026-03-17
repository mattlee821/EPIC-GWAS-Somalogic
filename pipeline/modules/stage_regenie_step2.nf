process STAGE_REGENIE_STEP2 {
  label        'low'
  publishDir   { "${params.outdir}/${study_id}/${protein_id}/${group}/008_regenie-step2" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), val(protein_id), val(chunk_id), val(chr), path(regenie_file)

  output:
  tuple val(study_id), val(group), val(protein_id), path("${chunk_id}_chr${chr}.regenie.gz")

  script:
  """
  cp ${regenie_file} ${chunk_id}_chr${chr}.regenie.gz
  """
}

process ADD_PCS {
  label        'low'
  conda        "${projectDir}/envs/r_gwas.yml"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/003_prepare-covariates" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), path(cov_file, stageAs: 'covariates_in.cov'), path(pcs_file), path(eval_file)

  output:
  tuple val(study_id), val(group), path("covariates.cov")
  path "scree_plot.png"

  script:
  """
  Rscript ${projectDir}/bin/merge_pcs.R \\
    --cov_file ${cov_file} \\
    --pcs_file ${pcs_file} \\
    --eval_file ${eval_file} \\
    --outdir .
  """
}

process VALIDATE_INPUTS {
  label        'low'
  publishDir   "${params.outdir}/_shared/001_validation", mode: 'copy'

  input:
  path samplesheet
  path phenotype_file, stageAs: 'pheno_input.txt'
  path covariate_file, stageAs: 'covar_input.txt'
  val include_studies
  val include_proteins
  val covariates

  output:
  path "validation_report.txt"

  script:
  """
  python ${projectDir}/bin/validate_inputs.py \\
    --samplesheet ${samplesheet} \\
    --phenotype_file pheno_input.txt \\
    --covariate_file covar_input.txt \\
    --include_studies "${include_studies}" \\
    --include_proteins "${include_proteins}" \\
    --covariates "${covariates}" \\
    --outdir .
  """
}

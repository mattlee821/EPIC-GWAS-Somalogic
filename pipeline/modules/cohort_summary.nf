process COHORT_SUMMARY {
  label        'low'
  conda        "${projectDir}/envs/r_gwas.yml"
  publishDir   "${params.outdir}", mode: 'copy'

  input:
  val metrics_files

  output:
  path "all_metrics.tsv"
  path("lambda_gc_distribution.png", optional: true)

  script:
  def metrics_arg = metrics_files.join(',')
  """
  Rscript ${projectDir}/bin/generate_cohort_summary.R \\
    --metrics_files "${metrics_arg}" \\
    --outdir .
  """
}

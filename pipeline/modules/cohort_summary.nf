process COHORT_SUMMARY {
  label        'low'
  conda        "${projectDir}/envs/r_gwas.yml"
  publishDir   "${params.outdir}", mode: 'copy'

  input:
  path metrics_list

  output:
  path "all_metrics.tsv"
  path("lambda_gc_distribution.png", optional: true)

  script:
  """
  Rscript ${projectDir}/bin/generate_cohort_summary.R \\
    --metrics_list ${metrics_list} \\
    --outdir .
  """

}

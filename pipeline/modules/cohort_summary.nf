process COHORT_SUMMARY {
  label        'low'
  conda        "${projectDir}/envs/py_validate.yml"
  publishDir   "${params.outdir}", mode: 'copy'

  input:
  path metrics_list

  output:
  path "all_metrics.tsv"

  script:
  """
  python ${projectDir}/bin/generate_cohort_summary.py \\
    --metrics_list ${metrics_list} \\
    --outdir .
  """

}

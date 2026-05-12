process SAMPLE_QC {
  label        'medium'
  conda        "${projectDir}/envs/plink2.yml"

  input:
  tuple val(study_id), val(pfile_prefix), val(mind), val(king_cutoff)

  output:
  tuple val(study_id), path("qc_pass_samples.txt")

  script:
  """
  bash ${projectDir}/bin/sample_qc.sh \\
    --pfile ${pfile_prefix} \\
    --study ${study_id} \\
    --mind ${mind} \\
    --king ${king_cutoff} \\
    --outdir .
  """
}

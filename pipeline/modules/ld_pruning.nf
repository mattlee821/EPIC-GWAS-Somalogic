process LD_PRUNING {
  label        'medium'
  conda        "${projectDir}/envs/plink2.yml"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/006_ld-pruning" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), val(bfile_prefix), path(qc_samples), path(qc_variants), val(ld_window_kb), val(ld_step), val(ld_r2)

  output:
  tuple val(study_id), val(group), path("step1_input.bed"), path("step1_input.bim"), path("step1_input.fam")

  script:
  """
  bash ${projectDir}/bin/ld_pruning.sh \\
    --bfile ${bfile_prefix} \\
    --keep ${qc_samples} \\
    --extract ${qc_variants} \\
    --study ${study_id} \\
    --window ${ld_window_kb} \\
    --step ${ld_step} \\
    --r2 ${ld_r2} \\
    --outdir .
  """
}

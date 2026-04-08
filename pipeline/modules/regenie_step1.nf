process REGENIE_STEP1 {
  label        'regenie_step1'
  conda        "${projectDir}/envs/regenie.yml"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/007_regenie-step1" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), path(step1_bed), path(step1_bim), path(step1_fam), path(pheno_full), path(cov_file), val(bsize)

  output:
  tuple val(study_id), val(group), path("pred.list"), path("loco_*.loco.gz")

  script:
  // Step 1 expects the bfile prefix
  def bfile_prefix = step1_bed.baseName
  """
  bash ${projectDir}/bin/run_regenie_step1.sh \\
    --step1_bfile ${bfile_prefix} \\
    --pheno_file ${pheno_full} \\
    --cov_file ${cov_file} \\
    --study ${study_id} \\
    --group ${group} \\
    --bsize ${bsize} \\
    --outdir .
  """

}

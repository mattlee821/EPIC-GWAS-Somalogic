process VARIANT_QC {
  label        'medium'
  conda        "${projectDir}/envs/plink2.yml"

  input:
  tuple val(study_id), val(bfile_prefix), path(qc_samples), val(maf), val(hwe), val(geno)

  output:
  tuple val(study_id), path("variant_qc_pass.snplist")

  script:
  """
  bash ${projectDir}/bin/variant_qc.sh \\
    --bfile ${bfile_prefix} \\
    --keep ${qc_samples} \\
    --study ${study_id} \\
    --maf ${maf} \\
    --hwe ${hwe} \\
    --geno ${geno} \\
    --outdir .
  """
}

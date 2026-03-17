process PCA {
  label        'medium'
  conda        "${projectDir}/envs/plink2.yml"

  input:
  tuple val(study_id), val(group), path(step1_bed), path(step1_bim), path(step1_fam)

  output:
  tuple val(study_id), val(group), path("pca.eigenvec"), path("pca.eigenval")

  script:
  def prefix = step1_bed.baseName
  """
  bash ${projectDir}/bin/compute_pcs.sh \\
    --bfile ${prefix} \\
    --outdir .
  """
}

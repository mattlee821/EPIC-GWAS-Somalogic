process REGENIE_STEP2 {
  label        'regenie_step2'
  conda        "${projectDir}/envs/regenie.yml"
  publishDir   params.outdir, mode: 'link', saveAs: { name ->
    def m = (name =~ /.*_chr(\d+|X)_.+\.regenie\.gz$/)
    def chr_id = m.size() > 0 ? m[0][1] : "UNKNOWN"
    "${study_id}/${protein_id}/${group}/008_regenie-step2/chr${chr_id}.regenie.gz"
  }

  input:
  tuple val(study_id), val(group), val(protein_id), val(plink_bfile), val(bgen_dir), val(sample_file), val(chromosomes), path(pheno_file), path(cov_file), path(pred_list), path(loco_files), val(bsize)

  output:
  tuple val(study_id), val(group), val(protein_id), path("*.regenie.gz")

  script:
  """
  bash ${projectDir}/bin/run_regenie_step2.sh \\
    --pfile ${plink_bfile} \\
    --bgen_dir ${bgen_dir} \\
    --sample_file ${sample_file} \\
    --pheno_file ${pheno_file} \\
    --cov_file ${cov_file} \\
    --pred_list ${pred_list} \\
    --study ${study_id} \\
    --group ${group} \\
    --chromosomes "${chromosomes}" \\
    --bsize ${bsize} \\
    --outdir .
  """
}

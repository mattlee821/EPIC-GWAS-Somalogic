process REGENIE_STEP2 {
  label        'regenie_step2'
  conda        "${projectDir}/envs/regenie.yml"
  publishDir   params.outdir, mode: 'link', saveAs: { name ->
    def m = (name =~ /.*_chr(\d+|X)_(.*)\.regenie.gz$/)
    def protein_id = m.size() > 0 ? m[0][2] : "UNKNOWN"
    "${study_id}/${protein_id}/${group}/008_regenie-step2/chr${chr}.regenie.gz"
  }

  input:
  tuple val(study_id), val(group), val(chr), val(gen_input), val(gen_type), val(sample_file), path(pheno_file), path(cov_file), path(pred_list), path(loco_files), val(bsize)

  output:
  tuple val(study_id), val(group), val(chr), path("*.regenie.gz")

  script:
  """
  if [[ "${gen_type}" == "pfile" ]]; then
      gen_arg="--pfile ${gen_input}"
  else
      gen_arg="--bgen_file ${gen_input} --sample_file ${sample_file}"
  fi

  bash ${projectDir}/bin/run_regenie_step2.sh \\
    \${gen_arg} \\
    --pheno_file ${pheno_file} \\
    --cov_file ${cov_file} \\
    --pred_list ${pred_list} \\
    --study ${study_id} \\
    --group ${group} \\
    --chr ${chr} \\
    --bsize ${bsize} \\
    --outdir .
  """
}

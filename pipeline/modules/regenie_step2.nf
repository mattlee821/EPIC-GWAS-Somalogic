process REGENIE_STEP2 {
  label        'regenie_step2'
  conda        "${projectDir}/envs/regenie.yml"

  input:
  tuple val(study_id), val(group), val(chunk_id), val(chr), val(gen_input), val(gen_type), val(sample_file), path(pheno_chunk), path(cov_file), path(pred_list), path(loco_files), val(bsize)

  output:
  tuple val(study_id), val(group), val(chunk_id), val(chr), path("*.regenie.gz")

  script:
  """
  if [[ "${gen_type}" == "pfile" ]]; then
      gen_arg="--pfile ${gen_input}"
  else
      gen_arg="--bgen_file ${gen_input} --sample_file ${sample_file}"
  fi

  bash ${projectDir}/bin/run_regenie_step2.sh \\
    \${gen_arg} \\
    --pheno_file ${pheno_chunk} \\
    --cov_file ${cov_file} \\
    --pred_list ${pred_list} \\
    --study ${study_id} \\
    --group ${group} \\
    --chunk_id ${chunk_id} \\
    --chr ${chr} \\
    --bsize ${bsize} \\
    --outdir .
  """
}

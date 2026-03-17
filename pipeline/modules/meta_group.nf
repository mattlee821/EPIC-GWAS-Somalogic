process META_GROUP {
  label        'medium'
  // uses system-installed METAL binary (tools/metal/metal)
  publishDir   { "${params.outdir}/${study_id}/${protein_id}/meta/METAL" }, mode: 'copy'

  input:
  tuple val(study_id), val(protein_id), val(sumstats_files)

  output:
  tuple val(study_id), val(protein_id), path("meta.tbl.gz"), path("meta.tbl.info"), path("meta.log")

  script:
  def files_arg = sumstats_files.join(',')
  def metal_bin = file(params.metal_binary).toAbsolutePath().toString()
  """
  bash ${projectDir}/bin/run_metal.sh \\
    --input_files "${files_arg}" \\
    --protein_id ${protein_id} \\
    --group meta \\
    --metal_binary "${metal_bin}" \\
    --outdir .
  """
}

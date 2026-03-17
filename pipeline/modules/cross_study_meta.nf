process CROSS_STUDY_META {
  label        'medium'
  // uses system-installed METAL binary (tools/metal/metal)
  publishDir   { "${params.outdir}/meta/${protein_id}/${group}/METAL" }, mode: 'copy'

  input:
  tuple val(group), val(protein_id), val(sumstats_files)

  output:
  tuple val(group), val(protein_id), path("meta.tbl.gz"), path("meta.tbl.info"), path("meta.log")

  script:
  def files_arg = sumstats_files.join(',')
  def metal_bin = file(params.metal_binary).toAbsolutePath().toString()
  """
  bash ${projectDir}/bin/run_metal.sh \\
    --input_files "${files_arg}" \\
    --protein_id ${protein_id} \\
    --group ${group} \\
    --metal_binary "${metal_bin}" \\
    --outdir .
  """
}

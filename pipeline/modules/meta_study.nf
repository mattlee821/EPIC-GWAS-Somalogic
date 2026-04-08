process META_STUDY {
  label        'medium'
  // uses system-installed METAL binary (tools/metal/metal)
  publishDir   {
    def pid = (protein_ids instanceof java.util.Collection && protein_ids) ? protein_ids[0] : protein_ids
    "${params.outdir}/meta/${pid}/${group}/012_METAL"
  }, mode: 'copy'

  input:
  tuple val(group), val(protein_ids), val(sumstats_files), val(mapping)

  output:
  tuple val("meta_study"), val(group), val(protein_ids), path("meta.tbl.gz"), path("meta.tbl.info"), path("meta.log"), emit: all_out
  path "meta.tbl.gz", emit: meta_files

  script:
  def metal_bin = file(params.metal_binary).toAbsolutePath().toString()
  """
  echo -e "${mapping}" > mapping.txt
  while IFS=':' read -r prot files; do
    bash ${projectDir}/bin/run_metal.sh \\
      --input_files "\${files}" \\
      --protein_id "\${prot}" \\
      --group "${group}" \\
      --metal_binary "${metal_bin}" \\
      --outdir .
  done < mapping.txt
  """

}

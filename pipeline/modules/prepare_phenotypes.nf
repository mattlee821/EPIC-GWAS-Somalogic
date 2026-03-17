process PREPARE_PHENOTYPES {
  label        'low'
  conda        "${projectDir}/envs/r_gwas.yml"
  publishDir   { "${params.outdir}/${study_id}/_shared/${group}/002_prepare-phenotypes" }, mode: 'copy'

  input:
  tuple val(study_id), val(group), path(phenotype_file, stageAs: 'phenotype_input.txt'), path(sample_file), path(covariate_file, stageAs: 'covariate_input.txt'), val(group_column), val(cases_value), val(include_proteins), val(chunk_size), val(covariates)

  output:
  tuple val(study_id), val(group), path("chunk_*.pheno", optional: true), path("full.pheno"), path("chunks.manifest"), path("keep_samples.txt")
  path "protein_summary.tsv"

  script:
  """
  Rscript ${projectDir}/bin/prepare_phenotypes.R \\
    --phenotype_file phenotype_input.txt \\
    --sample_file ${sample_file} \\
    --study_id ${study_id} \\
    --group "${group}" \\
    --group_column "${group_column}" \\
    --cases_value "${cases_value}" \\
    --covariate_file covariate_input.txt \\
    --covariates "${covariates}" \\
    --include_proteins "${include_proteins}" \\
    --chunk_size ${chunk_size} \\
    --outdir .
  """
}

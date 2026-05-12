#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { VALIDATE_INPUTS } from './modules/validate_inputs.nf'
include { PREPARE_PHENOTYPES } from './modules/prepare_phenotypes.nf'
include { PREPARE_COVARIATES } from './modules/prepare_covariates.nf'
include { SAMPLE_QC } from './modules/sample_qc.nf'
include { VARIANT_QC } from './modules/variant_qc.nf'
include { STAGE_SAMPLE_QC } from './modules/stage_sample_qc.nf'
include { STAGE_VARIANT_QC } from './modules/stage_variant_qc.nf'
include { LD_PRUNING } from './modules/ld_pruning.nf'
include { PCA } from './modules/pca.nf'
include { ADD_PCS } from './modules/add_pcs.nf'
include { REGENIE_STEP0 } from './modules/regenie_step0.nf'
include { REGENIE_STEP1 } from './modules/regenie_step1.nf'
include { REGENIE_STEP2 } from './modules/regenie_step2.nf'
include { GWAS_QC } from './modules/gwas_qc.nf'
include { META_STUDY } from './modules/meta_study.nf'
include { META_GROUP } from './modules/meta_group.nf'
include { COHORT_SUMMARY } from './modules/cohort_summary.nf'

workflow {

    // 1. Validation
    // Check if covariate file exists, otherwise fallback to phenotype file for staging
    def cov_file = params.covariate_file && file(params.covariate_file).exists() ? file(params.covariate_file) : file(params.phenotype_file)

    VALIDATE_INPUTS(
        file(params.samplesheet),
        file(params.phenotype_file),
        cov_file,
        params.include_studies,
        params.include_proteins,
        params.covariates,
        params.chromosomes.collect { it.toString() }.join(',')
    )

    // 2. Load samplesheet and create study channel
    def required_studies = params.include_studies ? params.include_studies.split(',').collect { it.trim() } : []

    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .filter { row -> required_studies.size() == 0 || required_studies.contains(row.study_id) }
        .set { ch_studies }

    def ch_selected_study_count = ch_studies
        .map { row -> row.study_id.toString() }
        .unique()
        .collect()
        .map { study_ids -> study_ids.size() }

    // 3. Prepare covariates and phenotypes (branched by analysis groups)
    def groups = params.analysis_groups
    def ch_groups = Channel.fromList(groups)

    ch_studies
        .combine(Channel.fromList(groups))
        .filter { row, group ->
            if (group == "combined") return true
            return (row.group_column && row.group_column.trim() != "")
        }
        .map { row, group ->
            tuple(
                row.study_id,
                group,
                file(params.phenotype_file),
                file(row.sample_file),
                cov_file,
                row.group_column,
                row.cases_value,
                params.include_proteins,
                params.covariates
            )
        }
        .set { ch_pheno_inputs }

    ch_studies
        .combine(Channel.fromList(groups))
        .filter { row, group ->
            if (group == "combined") return true
            return (row.group_column && row.group_column.trim() != "")
        }
        .map { row, group ->
            tuple(row.study_id, group, cov_file, file(row.sample_file), row.group_column, row.cases_value, params.covariates, params.include_proteins)
        }
        .set { ch_cov_inputs }

    PREPARE_PHENOTYPES(ch_pheno_inputs)
    PREPARE_COVARIATES(ch_cov_inputs)

    // 4. Genetic QC per study (independent of analysis groups)
    ch_studies.map { row -> tuple(row.study_id, row.pfile, params.mind, params.king_cutoff) } | SAMPLE_QC

    // combine SAMPLE_QC (sid, qc_pass_samples.txt) with original pfile
    ch_studies_with_qc = ch_studies.map { row -> tuple(row.study_id, row.pfile) }
                                   .join(SAMPLE_QC.out, by: 0) // [sid, pfile, pass_samples]

    // VARIANT_QC takes [sid, pfile, qc_samples, maf, hwe, geno]
    VARIANT_QC(ch_studies_with_qc.map { sid, pfile, qcs -> tuple(sid, pfile, qcs, params.maf, params.hwe, params.geno) })

    // Stage study-level QC outputs into each group directory
    STAGE_SAMPLE_QC(
        SAMPLE_QC.out.combine(ch_groups)
                    .map { sid, qc_file, grp -> tuple(sid, grp, qc_file) }
    )
    STAGE_VARIANT_QC(
        VARIANT_QC.out.combine(ch_groups)
                     .map { sid, qc_var, grp -> tuple(sid, grp, qc_var) }
    )

    // 4b. Identify intersection of QCed samples and phenotyped samples for LD pruning
    pheno_samples = PREPARE_PHENOTYPES.out[0].map { sid, grp, feature_files, full, manifest, keep -> tuple(sid, grp, keep) }

    ch_for_ld = ch_studies_with_qc.join(VARIANT_QC.out, by: 0)
                                  .combine(pheno_samples, by: 0)
                                  .map { sid, pfile, qcs_g, qcv_g, grp, qcs_p ->
                                      tuple(sid, grp, pfile, qcs_p, qcv_g, params.ld_window_kb, params.ld_step, params.ld_r2)
                                  }

    LD_PRUNING(ch_for_ld)

    // 4c. Compute PCs on the step1 input and merge into covariates
    PCA(LD_PRUNING.out.map { sid, grp, pgen, pvar, psam -> tuple(sid, grp, pgen, pvar, psam) })

    def cov_for_pcs = PREPARE_COVARIATES.out.map { sid, grp, cov -> tuple("${sid}_${grp}", sid, grp, cov) }
    def pcs_for_cov = PCA.out.map { sid, grp, pvec, pval -> tuple("${sid}_${grp}", pvec, pval) }

    def cov_with_pcs_in = cov_for_pcs.join(pcs_for_cov, by: 0)
                                     .map { key, sid, grp, cov, pvec, pval -> tuple(sid, grp, cov, pvec, pval) }

    ADD_PCS(cov_with_pcs_in)

    // 5. Expand to one phenotype file per study, group, and feature
    def pheno_features = PREPARE_PHENOTYPES.out[0]
        .map { sid, grp, feature_files, full, manifest, keep -> tuple(sid, grp, feature_files, manifest) }
        .splitCsv(elem: 3, header: true, sep: '\t')
        .map { row ->
            def protein_id = row[3].protein_id
            def pheno_name = new File(row[3].pheno_file).name
            def feature_files = row[2] instanceof java.util.Collection ? row[2] : [row[2]]
            def pheno_file = feature_files.find { it.name == pheno_name }
            if (pheno_file == null) {
                throw new IllegalStateException("Could not resolve phenotype file '${pheno_name}' for ${row[0]} ${row[1]} ${protein_id}.")
            }
            tuple(row[0], row[1], protein_id, pheno_file)
        }

    // 6. REGENIE Step 0 (shared level-0 work per study and group)
    def regenie_key = { sid, grp -> "${sid}__${grp}" }
    def regenie_full = PREPARE_PHENOTYPES.out[0].map { sid, grp, feature_files, full, manifest, keep -> tuple(regenie_key(sid, grp), full) }
    def regenie_cov = ADD_PCS.out[0].map { sid, grp, cov -> tuple(regenie_key(sid, grp), cov) }
    def regenie_ld = LD_PRUNING.out.map { sid, grp, pgen, pvar, psam -> tuple(regenie_key(sid, grp), sid, grp, pgen, pvar, psam) }

    def regenie_shared = regenie_ld
        .join(regenie_full, by: 0)
        .join(regenie_cov, by: 0)
        .map { key, sid, grp, pgen, pvar, psam, full, cov ->
            tuple(key, sid, grp, pgen, pvar, psam, full, cov)
        }

    step0_inputs = regenie_shared
        .map { key, sid, grp, pgen, pvar, psam, full, cov ->
            tuple(sid, grp, pgen, pvar, psam, full, cov, params.regenie_bsize_step1)
        }

    REGENIE_STEP0(step0_inputs)

    // 7. REGENIE Step 1 (level-1 models per study, group, and feature)
    def pheno_prep = pheno_features.map { sid, grp, protein_id, pheno_file -> tuple(regenie_key(sid, grp), sid, grp, protein_id) }
    def step0_master = REGENIE_STEP0.out.map { sid, grp, master, snplists, l0_files ->
        tuple(regenie_key(sid, grp), master.toString())
    }

    step1_inputs = pheno_prep
        .combine(regenie_shared.map { key, sid, grp, pgen, pvar, psam, full, cov -> tuple(key, pgen, pvar, psam, full, cov) }, by: 0)
        .combine(step0_master, by: 0)
        .map { key, sid, grp, protein_id, pgen, pvar, psam, full, cov, master_file ->
            tuple(sid, grp, protein_id, pgen, pvar, psam, full, cov, master_file, params.regenie_bsize_step1)
        }

    REGENIE_STEP1(step1_inputs)

    // 8. REGENIE Step 2 (one run per study, group, and feature across all chromosomes)
    def chromosome_list = params.chromosomes.collect { it.toString() }.join(',')
    def ch_step1 = REGENIE_STEP1.out.map { sid, grp, protein_id, pred, loco -> tuple("${sid}__${grp}__${protein_id}", pred, loco) }
    def ch_pheno_step2 = pheno_features.map { sid, grp, protein_id, pheno_file -> tuple("${sid}__${grp}__${protein_id}", sid, grp, protein_id, pheno_file) }
    def ch_cov_step2 = ADD_PCS.out[0].map { sid, grp, cov -> tuple("${sid}_${grp}", cov) }
    def ch_genotype = ch_studies.map { r -> tuple(r.study_id, r.pfile, r.sample_file) }

    step2_setup = ch_pheno_step2
        .join(ch_step1, by: 0)
        .map { key, sid, grp, protein_id, pheno_file, pred, loco -> tuple("${sid}_${grp}", sid, grp, protein_id, pheno_file, pred, loco) }
        .combine(ch_cov_step2, by: 0)
        .map { key, sid, grp, protein_id, pheno_file, pred, loco, cov -> tuple(sid, grp, protein_id, pheno_file, pred, loco, cov) }
        .combine(ch_genotype, by: 0)
        .map { sid, grp, protein_id, pheno_file, pred, loco, cov, pfile, samp ->
            tuple(sid, grp, protein_id, pfile, samp, chromosome_list, pheno_file, cov, pred, loco, params.regenie_bsize_step2)
        }

    REGENIE_STEP2(step2_setup)

    // 9. Group REGENIE step 2 outputs by protein
    ch_grouped_regenie = REGENIE_STEP2.out
        .map { sid, grp, protein_id, files ->
            def result_files = files instanceof java.util.Collection ? files : [files]
            tuple(sid, grp, protein_id, result_files)
        }

    // 10. Consolidated post-GWAS QC
    def gwas_qc_batch_size = (params.gwas_qc_batch_size ?: 50) as int
    if (gwas_qc_batch_size < 1) {
        throw new IllegalArgumentException("params.gwas_qc_batch_size must be >= 1")
    }

    ch_gwas_qc_in = ch_grouped_regenie
        .map { sid, grp, prot, fs ->
            def files = (fs instanceof java.util.Collection ? fs : [fs])
            tuple("${sid}__${grp}", sid, grp, prot.toString(), files)
        }
        .groupTuple(by: 0, size: gwas_qc_batch_size, remainder: true)
        .map { key, sids, grps, prots, file_sets ->
            def prot_ids = prots.collect { it.toString() }
            def files = file_sets.collectMany { it instanceof java.util.Collection ? it : [it] }
            def mapping = (0..<prot_ids.size()).collect { idx ->
                def prot_files = file_sets[idx] instanceof java.util.Collection ? file_sets[idx] : [file_sets[idx]]
                "${prot_ids[idx]}:${prot_files.collect { it.name }.join(',')}"
            }.join('\n')
            tuple(sids[0].toString(), grps[0].toString(), prot_ids, files, mapping)
        }

    GWAS_QC(ch_gwas_qc_in)

    // Flatten batched GWAS_QC outputs for downstream use
    ch_gwas_qc_results = GWAS_QC.out.all_out
        .flatMap { sid, grp, prots, gwas_files, metrics_files ->
            def prot_ids = prots instanceof java.util.Collection ? prots.collect { it.toString() } : [prots.toString()]
            def gwas_list = gwas_files instanceof java.util.Collection ? gwas_files : [gwas_files]
            def metrics_list = metrics_files instanceof java.util.Collection ? metrics_files : [metrics_files]

            def protein_from_output = { f ->
                def parent = f.getParent()?.getFileName()?.toString()
                if (parent) {
                    return parent
                }
                def parts = f.toString().tokenize('/\\')
                return parts.size() >= 2 ? parts[-2] : null
            }

            def gwas_by_prot = gwas_list.collectEntries { f ->
                [(protein_from_output(f)): f]
            }
            def metrics_by_prot = metrics_list.collectEntries { f ->
                [(protein_from_output(f)): f]
            }

            prot_ids.collect { prot_id ->
                def gwas = gwas_by_prot[prot_id]
                def metrics = metrics_by_prot[prot_id]
                if (gwas == null || metrics == null) {
                    throw new IllegalStateException(
                        "Could not resolve GWAS_QC outputs for ${sid} ${grp} ${prot_id}. " +
                        "GWAS keys=${gwas_by_prot.keySet().findAll { it != null }.sort().join(',')} " +
                        "Metrics keys=${metrics_by_prot.keySet().findAll { it != null }.sort().join(',')}"
                    )
                }
                tuple(sid, grp, prot_id, gwas, metrics)
            }
        }

    // 8b. Meta-analysis of cases vs controls within each study (optional)
    def meta_group_enabled = params.meta_group?.toString()?.toLowerCase() in ['true', '1', 'yes']
    def ch_meta_group_qc_out = Channel.empty()
    def ch_meta_group_for_meta = Channel.empty()

    if (meta_group_enabled) {
        def ch_meta_group_in = ch_gwas_qc_results
            .filter { sid, grp, prot, gwas, metrics -> grp == "cases" || grp == "controls" }
            .map { sid, grp, prot, gwas, metrics -> tuple("${sid}__${prot}", sid, prot, gwas) }
            .groupTuple(by: 0, size: 2, remainder: true)
            .filter { key, sids, prots, gwases -> gwases.size() == 2 }
            .map { key, sids, prots, gwases ->
                def files = gwases instanceof java.util.Collection ? gwases : [gwases]
                tuple(sids[0].toString(), prots[0].toString(), files)
            }

        META_GROUP(ch_meta_group_in)

        ch_meta_group_qc_out = META_GROUP.out.all_out
        ch_meta_group_for_meta = META_GROUP.out.all_out
            .map { sid, grp, prot, gwas, met -> tuple(sid, grp, prot, gwas) }
    }

    // 9. Meta-analysis across studies
    def meta_study_raw = params.meta_study
    def meta_study_enabled = meta_study_raw && !(meta_study_raw.toString().trim().toLowerCase() in ['false', '0', 'no'])

    if (meta_study_enabled) {
        def meta_study_list = meta_study_raw.toString().split(',').collect { it.trim().toLowerCase() }.findAll { it }
        if (meta_study_list == ["true"] || meta_study_list == ["1"]) meta_study_list = ["meta"]
        def meta_study_set = meta_study_list as Set

        ch_gwas_qc_results.map { sid, grp, prot, gwas, met -> tuple(sid, grp, prot, gwas) }
            .mix(ch_meta_group_for_meta)
            .filter { sid, grp, prot, gwas -> meta_study_set.contains(grp.toString().toLowerCase()) }
            .combine(ch_selected_study_count)
            .map { sid, grp, prot, gwas, study_count ->
                tuple(groupKey("${grp}__${prot}", study_count as int), sid, grp, prot, gwas)
            }
            .groupTuple(remainder: true)
            .filter { key, sids, grps, prots, gwases -> gwases.size() > 1 }
            .map { key, sids, grps, prots, gwases ->
                def files = gwases instanceof java.util.Collection ? gwases : [gwases]
                tuple(grps[0].toString(), prots[0].toString(), files)
            }
            .set { ch_meta_study_in }

        META_STUDY(ch_meta_study_in)
    }

    // 10. Cohort summary
    ch_metrics_all = ch_gwas_qc_results.map { sid, grp, prot, gwas, met -> met.toString() }

    if (meta_group_enabled) {
        ch_metrics_all = ch_metrics_all.mix(ch_meta_group_qc_out.map { sid, grp, prot, gwas, met -> met.toString() })
    }

    if (meta_study_enabled) {
        ch_metrics_all = ch_metrics_all.mix(META_STUDY.out.all_out.map { sid, grp, prot, gwas, met -> met.toString() })
    }

    def ch_metrics_list = ch_metrics_all.collectFile(name: 'metrics_list.txt', newLine: true)
    COHORT_SUMMARY(ch_metrics_list)
}

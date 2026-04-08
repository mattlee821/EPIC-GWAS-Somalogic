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
include { REGENIE_STEP1 } from './modules/regenie_step1.nf'
include { REGENIE_STEP2 } from './modules/regenie_step2.nf'
include { GWAS_QC } from './modules/gwas_qc.nf'
include { META_STUDY } from './modules/meta_study.nf'
include { META_GROUP } from './modules/meta_group.nf'
include { META_QC as META_QC_GROUP } from './modules/meta_qc.nf'
include { META_QC as META_QC_STUDY } from './modules/meta_qc.nf'
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
        params.covariates
    )

    // 2. Load samplesheet and create study channel
    def required_studies = params.include_studies ? params.include_studies.split(',').collect { it.trim() } : []

    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .filter { row -> required_studies.size() == 0 || required_studies.contains(row.study_id) }
        .set { ch_studies }

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
    ch_studies.map { row -> tuple(row.study_id, row.plink_bfile, params.mind, params.king_cutoff) } | SAMPLE_QC

    // combine SAMPLE_QC (sid, qc_pass_samples.txt) with original bfile
    ch_studies_with_qc = ch_studies.map { row -> tuple(row.study_id, row.plink_bfile) }
                                   .join(SAMPLE_QC.out, by: 0) // [sid, bfile, pass_samples]

    // VARIANT_QC takes [sid, bfile, qc_samples, maf, hwe, geno]
    VARIANT_QC(ch_studies_with_qc.map { sid, bfile, qcs -> tuple(sid, bfile, qcs, params.maf, params.hwe, params.geno) })

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
    pheno_samples = PREPARE_PHENOTYPES.out[0].map { sid, grp, full, keep -> tuple(sid, grp, keep) }

    ch_for_ld = ch_studies_with_qc.join(VARIANT_QC.out, by: 0)
                                  .combine(pheno_samples, by: 0)
                                  .map { sid, pfile, qcs_g, qcv_g, grp, qcs_p ->
                                      tuple(sid, grp, pfile, qcs_p, qcv_g, params.ld_window_kb, params.ld_step, params.ld_r2)
                                  }

    LD_PRUNING(ch_for_ld)

    // 4c. Compute PCs on the step1 input and merge into covariates
    PCA(LD_PRUNING.out.map { sid, grp, bed, bim, fam -> tuple(sid, grp, bed, bim, fam) })

    def cov_for_pcs = PREPARE_COVARIATES.out.map { sid, grp, cov -> tuple("${sid}_${grp}", sid, grp, cov) }
    def pcs_for_cov = PCA.out.map { sid, grp, pvec, pval -> tuple("${sid}_${grp}", pvec, pval) }

    def cov_with_pcs_in = cov_for_pcs.join(pcs_for_cov, by: 0)
                                     .map { key, sid, grp, cov, pvec, pval -> tuple(sid, grp, cov, pvec, pval) }

    ADD_PCS(cov_with_pcs_in)

    // 5. REGENIE Step 1 (one run per study and group)
    def pheno_prep = PREPARE_PHENOTYPES.out[0].map { sid, grp, full, keep -> tuple("${sid}_${grp}", full) }
    def cov_prep = ADD_PCS.out[0].map { sid, grp, cov -> tuple("${sid}_${grp}", cov) }
    def ld_prep = LD_PRUNING.out.map { sid, grp, bed, bim, fam -> tuple("${sid}_${grp}", sid, grp, bed, bim, fam) }

    step1_inputs = ld_prep.join(pheno_prep, by: 0)
                          .join(cov_prep, by: 0)
                          .map { key, sid, grp, bed, bim, fam, full, cov ->
                              tuple(sid, grp, bed, bim, fam, full, cov, params.regenie_bsize_step1)
                          }

    REGENIE_STEP1(step1_inputs)

    // 6. REGENIE Step 2 (one run per study, group, chromosome)
    def chrs_ch = Channel.fromList(params.chromosomes)
    def ch_step1 = REGENIE_STEP1.out.map { sid, grp, pred, loco -> tuple(sid, grp, pred, loco) }
    def ch_pheno_step2 = PREPARE_PHENOTYPES.out[0].map { sid, grp, full, keep -> tuple(sid, grp, full) }
    def ch_cov_step2 = ADD_PCS.out[0].map { sid, grp, cov -> tuple(sid, grp, cov) }
    def ch_genotype = ch_studies.map { r -> tuple(r.study_id, r.plink_bfile, r.bgen_dir, r.sample_file) }

    step2_setup = ch_step1
        .join(ch_pheno_step2, by: [0, 1])
        .join(ch_cov_step2, by: [0, 1])
        .combine(chrs_ch)
        .map { sid, grp, pred, loco, full, cov, chr -> tuple(sid, grp, chr, pred, loco, full, cov) }
        .combine(ch_genotype, by: 0)
        .map { sid, grp, chr, pred, loco, full, cov, pbfile, bgendir, samp ->
            def bgen_file = file("${bgendir}/chr${chr}.bgen")
            def gen_input
            def gen_type
            if (bgen_file.exists()) {
                gen_input = bgen_file.toString()
                gen_type = "bgen"
            } else {
                gen_input = pbfile
                gen_type = "pfile"
            }
            tuple(sid, grp, chr, gen_input, gen_type, samp, full, cov, pred, loco, params.regenie_bsize_step2)
        }

    REGENIE_STEP2(step2_setup)

    // 7. Parse REGENIE step 2 outputs by protein
    ch_regenie_parsed = REGENIE_STEP2.out
        .flatMap { sid, grp, chr, files ->
            (files instanceof java.util.Collection ? files : [files]).collect { f ->
                def matcher = (f.name =~ /.*_chr(\d+|X)_(.*)\.regenie.gz$/)
                def prot = matcher.size() > 0 ? matcher[0][2] : "UNKNOWN"
                tuple(sid, grp, prot, chr, f)
            }
        }

    // Group REGENIE outputs by (study, group, protein)
    ch_grouped_regenie = ch_regenie_parsed
        .map { sid, grp, prot, chr_id, f -> tuple("${sid}__${grp}__${prot}", sid, grp, prot, f) }
        .groupTuple(by: 0)
        .map { key, sids, grps, prots, fs -> tuple(sids[0], grps[0], prots[0], fs) }

    // 8. Consolidated post-GWAS QC
    ch_gwas_qc_in = ch_grouped_regenie
        .map { sid, grp, prot, fs ->
            def files = (fs instanceof java.util.Collection ? fs : [fs])
            def mapping = "${prot}:${files.collect{ it.name }.join(',')}"
            tuple(sid, grp, prot, files, mapping)
        }

    GWAS_QC(ch_gwas_qc_in)

    // Flatten GWAS_QC outputs for downstream use
    ch_gwas_qc_results = GWAS_QC.out.all_out
        .transpose()
        .map { sid, grp, prot, gwas, metrics ->
            def prot_id = (prot instanceof java.util.Collection && prot) ? prot[0].toString() : prot.toString()
            tuple(sid, grp, prot_id, gwas, metrics)
        }

    // 8b. Meta-analysis of cases vs controls within each study (optional)
    def meta_group_enabled = params.meta_group?.toString()?.toLowerCase() in ['true', '1', 'yes']
    def ch_meta_group_qc_out = Channel.empty()
    def ch_meta_group_for_meta = Channel.empty()

    if (meta_group_enabled) {
        def ch_meta_group_prep = ch_gwas_qc_results
            .filter { sid, grp, prot, gwas, metrics -> grp == "cases" || grp == "controls" }
            .map { sid, grp, prot, gwas, metrics -> tuple("${sid}__${prot}", sid, prot, gwas) }
            .groupTuple(by: 0)
            .filter { key, sids, prots, gwases -> gwases.size() == 2 }
            .map { key, sids, prots, gwases -> tuple(sids[0], prots[0], gwases) }

        def ch_meta_group_in = ch_meta_group_prep
            .map { sid, prot, gwases ->
                def files = (gwases instanceof java.util.Collection ? gwases : [gwases])
                def mapping = "${prot}:${files.collect { it.toString() }.join(',')}"
                tuple(sid, prot, files, mapping)
            }

        META_GROUP(ch_meta_group_in)

        def ch_meta_group_qc_in = META_GROUP.out.all_out
            .transpose()
            .map { sid, prot, meta_tbl, meta_info, meta_log ->
                def prot_id = (prot instanceof java.util.Collection && prot) ? prot[0].toString() : prot.toString()
                tuple(sid, "meta", prot_id, meta_tbl, "${prot_id}:${meta_tbl.name}")
            }

        META_QC_GROUP(ch_meta_group_qc_in)

        ch_meta_group_qc_out = META_QC_GROUP.out.all_out.transpose()
        ch_meta_group_for_meta = ch_meta_group_qc_out.map { sid, grp, prot, gwas, met -> tuple(sid, grp, prot, gwas) }
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
            .map { sid, grp, prot, gwas -> tuple("${grp}__${prot}", sid, grp, prot, gwas) }
            .groupTuple(by: 0)
            .filter { key, sids, grps, prots, gwases -> gwases.size() > 1 }
            .map { key, sids, grps, prots, gwases -> tuple(grps[0], prots[0], gwases) }
            .map { grp, prot, gwases ->
                def files = (gwases instanceof java.util.Collection ? gwases : [gwases])
                def mapping = "${prot}:${files.collect { it.toString() }.join(',')}"
                tuple(grp, prot, files, mapping)
            }
            .set { ch_meta_study_in }

        META_STUDY(ch_meta_study_in)

        ch_meta_study_qc_in = META_STUDY.out.all_out
            .transpose()
            .map { sid, grp, prot, meta_tbl, meta_info, meta_log ->
                def prot_id = (prot instanceof java.util.Collection && prot) ? prot[0].toString() : prot.toString()
                tuple("meta_study", grp, prot_id, meta_tbl, "${prot_id}:${meta_tbl.name}")
            }

        META_QC_STUDY(ch_meta_study_qc_in)
    }

    // 10. Cohort summary
    ch_metrics_all = ch_gwas_qc_results.map { sid, grp, prot, gwas, met -> met.toString() }

    if (meta_group_enabled) {
        ch_metrics_all = ch_metrics_all.mix(ch_meta_group_qc_out.map { sid, grp, prot, gwas, met -> met.toString() })
    }

    if (meta_study_enabled) {
        def ch_meta_study_qc_all = META_QC_STUDY.out.all_out.transpose()
        ch_metrics_all = ch_metrics_all.mix(ch_meta_study_qc_all.map { sid, grp, prot, gwas, met -> met.toString() })
    }

    def ch_metrics_list = ch_metrics_all.collectFile(name: 'metrics_list.txt', newLine: true)
    COHORT_SUMMARY(ch_metrics_list)
}

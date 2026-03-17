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
include { STAGE_REGENIE_STEP2 } from './modules/stage_regenie_step2.nf'
include { POST_GWAS_QC } from './modules/post_gwas_qc.nf'
include { CROSS_STUDY_META } from './modules/cross_study_meta.nf'
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
    def required_studies = params.include_studies ? params.include_studies.split(',').collect{it.trim()} : []
    
    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .filter { row -> required_studies.size() == 0 || required_studies.contains(row.study_id) }
        .set { ch_studies_raw }

    // Count studies to decide on meta-analysis
    ch_studies_count = ch_studies_raw.count()
    ch_studies = ch_studies_raw
        
    // 3. Prepare Covariates and Phenotypes (branched by analysis groups)
    def groups = params.analysis_groups
    def ch_groups = Channel.fromList(groups)
    
    ch_studies
        .combine(Channel.fromList(groups))
        .filter { row, group -> 
            if (group == "combined") return true
            return (row.group_column && row.group_column.trim() != "")
        }
        .map { row, group -> 
            tuple(row.study_id, group, file(params.phenotype_file), file(row.sample_file), cov_file, row.group_column, row.cases_value, params.include_proteins, params.chunk_size, params.covariates)
        }.set { ch_pheno_inputs }
        
    ch_studies
        .combine(Channel.fromList(groups))
        .filter { row, group -> 
            if (group == "combined") return true
            return (row.group_column && row.group_column.trim() != "")
        }
        .map { row, group ->
            tuple(row.study_id, group, cov_file, file(row.sample_file), row.group_column, row.cases_value, params.include_covariates, params.include_proteins)
        }.set { ch_cov_inputs }

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
    
    // 4.  Identify intersection of QCed samples and phenotyped samples for LD pruning
    // PREPARE_PHENOTYPES out: tuple(sid, grp, chunks, full, man, keep_samples)
    // SAMPLE_QC out: tuple(sid, qc_pass_samples.txt)
    // VARIANT_QC out: tuple(sid, qc_variants.snplist)

    pheno_samples = PREPARE_PHENOTYPES.out[0].map { sid, grp, chunks, full, man, keep -> tuple(sid, grp, keep) }
    
    // Join study-level QC info with group-level pheno info
    ch_for_ld = ch_studies_with_qc.join(VARIANT_QC.out, by: 0)
                                     .combine(pheno_samples, by: 0) // tuple(sid, pfile, qcs_global, qcv_global, grp, qcs_pheno)
                                     .map { sid, pfile, qcs_g, qcv_g, grp, qcs_p -> 
                                         tuple(sid, grp, pfile, qcs_p, qcv_g, params.ld_window_kb, params.ld_step, params.ld_r2) 
                                     }

    LD_PRUNING(ch_for_ld)

    // 4b. Compute PCs on the step1 input and merge into covariates
    PCA(LD_PRUNING.out.map { sid, grp, bed, bim, fam -> tuple(sid, grp, bed, bim, fam) })

    def cov_for_pcs = PREPARE_COVARIATES.out.map { sid, grp, cov -> tuple("${sid}_${grp}", sid, grp, cov) }
    def pcs_for_cov = PCA.out.map { sid, grp, pvec, pval -> tuple("${sid}_${grp}", pvec, pval) }

    def cov_with_pcs_in = cov_for_pcs.join(pcs_for_cov, by: 0)
                                   .map { key, sid, grp, cov, pvec, pval -> tuple(sid, grp, cov, pvec, pval) }

    ADD_PCS(cov_with_pcs_in)
    
    // 5. REGENIE Step 1 (per study, per group)
    // Join phenotype/covariate info with its specific LD output
    // PREPARE_PHENOTYPES.out: [sid, grp, chunks, full, man, keep]
    // PREPARE_COVARIATES.out: [sid, grp, cov]
    // LD_PRUNING.out: [sid, grp, bed, bim, fam]

    def pheno_prep = PREPARE_PHENOTYPES.out[0].map { sid, grp, chunks, full, man, keep -> tuple("${sid}_${grp}", full) }
    def cov_prep = ADD_PCS.out[0].map { sid, grp, cov -> tuple("${sid}_${grp}", cov) }
    def ld_prep = LD_PRUNING.out.map { sid, grp, bed, bim, fam -> tuple("${sid}_${grp}", sid, grp, bed, bim, fam) }

    step1_inputs = ld_prep.join(pheno_prep, by: 0)
                          .join(cov_prep, by: 0)
                          .map { key, sid, grp, bed, bim, fam, full, cov -> 
                              tuple(sid, grp, bed, bim, fam, full, cov, params.regenie_bsize_step1) 
                          }
                            
    REGENIE_STEP1(step1_inputs)

    // 6. REGENIE Step 2 (per study, per group, per chunk, per chromosome)
    // man files contain list of chunks. We need to split them out into items.
    
    def chrs = Channel.fromList(params.chromosomes)
    
    // REGENIE_STEP1 out: [sid, grp, pred, loco]
    // REGENIE_STEP1 out: [sid, grp, pred, loco]
    def step2_prep = PREPARE_PHENOTYPES.out[0].map { sid, grp, chunks, full, man, keep -> tuple("${sid}_${grp}", man, chunks) }
    def step2_regenie = REGENIE_STEP1.out.map { sid, grp, pred, loco -> tuple("${sid}_${grp}", sid, grp, pred, loco) }

    step2_setup = step2_prep.join(step2_regenie, by: 0)
                            .map { k, man, chunks, sid, grp, pred, loco -> 
                                tuple(sid, tuple(grp, man, chunks, pred, loco)) 
                            }
                            .combine(ch_studies.map{ r -> tuple(r.study_id, r.plink_bfile, r.bgen_dir, r.sample_file) }, by: 0)
                            .map { sid, tup, pfile, bgen_dir, samp -> 
                                tuple(sid, tup[0], tup[1], tup[2], pfile, bgen_dir, samp, tup[3], tup[4]) 
                            }
                           // [sid, grp, man, chunks, pfile, bgen_dir, samp, pred, loco]
                           
    // To chunk, we parse the manifest. Nextflow splitCsv is perfect here.
    def chunk_manifest = step2_setup.splitCsv(elem: 2, header: true, sep: '\t')
                              .map { row -> 
                                  // row[2] is now the parsed dictionary from manifest
                                  def chunk_id = row[2].chunk_id
                                  def chunk_file_name = new File(row[2].pheno_file).name
                                  def prot_list = (row[2].proteins ?: "").split(',').findAll { it }
                                  // Find the actual file object from the list of chunks `row[3]`
                                  def f_list = (row[3] instanceof java.util.Collection) ? row[3] : [row[3]]
                                  def chnk_file = f_list.find { it.name == chunk_file_name }
                                  tuple(row[0], row[1], chunk_id, chnk_file, row[4], row[5], row[6], row[7], row[8], prot_list)
                              }

    step2_chunks = chunk_manifest
                              .combine(chrs)
                              .map { sid, grp, chk_id, chk_f, pfile, bgendir, samp, pred, loco, prot_list, chr -> 
                                  def bgen_file = file("${bgendir}/chr${chr}.bgen")
                                  def gen_input
                                  def gen_type
                                  
                                  if (bgen_file.exists()) {
                                      gen_input = bgen_file.toString()
                                      gen_type = "bgen"
                                  } else {
                                      // Fallback to pfile (prefix)
                                      gen_input = pfile
                                      gen_type = "pfile"
                                  }
                                  tuple(sid, grp, chk_id, chr, gen_input, gen_type, samp, chk_f, pred, loco)
                              }

    def chunk_prots = chunk_manifest
                              .map { sid, grp, chk_id, chk_f, pfile, bgendir, samp, pred, loco, prot_list -> tuple("${sid}__${grp}__${chk_id}", prot_list) }
                              .distinct()

    // Recover cov_file for step 2
    def cov_for_step2 = ADD_PCS.out[0].map{sid, grp, cov -> tuple("${sid}_${grp}", cov)}
    
    step2_final = step2_chunks.map{ sid, grp, chk_id, chr, gen, type, samp, chk_f, pred, loco -> tuple("${sid}_${grp}", sid, grp, chk_id, chr, gen, type, samp, chk_f, pred, loco) }
                              .combine(cov_for_step2, by: 0)
                              .map { k, sid, grp, chk_id, chr, gen, type, samp, chk_f, pred, loco, cov -> 
                                  tuple(sid, grp, chk_id, chr, gen, type, samp, chk_f, cov, pred, loco, params.regenie_bsize_step2) 
                              }
                              
    REGENIE_STEP2(step2_final)
    
    // 7. Parse REGENIE step 2 outputs by protein and stage raw files
    ch_regenie_parsed = REGENIE_STEP2.out
        .map { sid, grp, chk, chr, files -> tuple("${sid}__${grp}__${chk}", sid, grp, chk, chr, files) }
        .combine(chunk_prots, by: 0)
        .flatMap { key, sid, grp, chk, chr, files, prot_list -> 
            (files instanceof java.util.Collection ? files : [files]).collect { f -> 
                def matcher = (f.name =~ /.*_chr(\\d+|X)_(.*)\\.regenie.gz$/)
                def prot = matcher.size() > 0 ? matcher[0][2] : (prot_list.size() == 1 ? prot_list[0] : "UNKNOWN")
                tuple(sid, grp, prot, chk, chr, f)
            }
        }

    STAGE_REGENIE_STEP2(ch_regenie_parsed)

    // Group REGENIE outputs by (study, group, protein)
    ch_grouped_regenie = ch_regenie_parsed
        .map { sid, grp, prot, chunk_id, chr_id, f -> tuple("${sid}__${grp}__${prot}", sid, grp, prot, f) }
        .groupTuple(by: 0)
        .map { key, sids, grps, prots, fs -> tuple(sids[0], grps[0], prots[0], fs) }

    // 8. Consolidated post-GWAS QC (finalize, QC, plots)
    POST_GWAS_QC(ch_grouped_regenie)

    // 8b. Meta-analysis of cases vs controls within each study (optional)
    def meta_group_enabled = params.meta_group?.toString()?.toLowerCase() in ['true','1','yes']
    def ch_meta_group_for_meta = Channel.empty()
    def ch_meta_group_qc_out = Channel.empty()
    if (meta_group_enabled) {
        def ch_meta_group = POST_GWAS_QC.out
            .filter { sid, grp, prot, gwas, metrics, fig -> grp == "cases" || grp == "controls" }
            .map { sid, grp, prot, gwas, metrics, fig -> tuple("${sid}__${prot}", sid, prot, gwas.toString()) }
            .groupTuple(by: 0)
            .map { key, sids, prots, gwases -> tuple(sids[0], prots[0], gwases) }
            .filter { sid, prot, gwases -> gwases.size() == 2 }

        META_GROUP(ch_meta_group)
        def ch_meta_group_qc_in = META_GROUP.out.map { sid, prot, meta_tbl, meta_info, meta_log -> tuple(sid, "meta", prot, meta_tbl) }
        META_QC_GROUP(ch_meta_group_qc_in)
        ch_meta_group_qc_out = META_QC_GROUP.out
        ch_meta_group_for_meta = ch_meta_group_qc_out.map { sid, grp, prot, gwas, metrics, fig -> tuple(sid, grp, prot, gwas) }
    }
    
    // 9. Meta-analysis (Skip if only 1 study)
    // Meta-analysis across studies for selected analysis groups
    def meta_study_raw = params.meta_study
    def meta_study_enabled = meta_study_raw && !(meta_study_raw.toString().trim().toLowerCase() in ['false','0','no'])
    log.info "Meta-study selection: ${meta_study_raw ?: 'false'}"
    if (!meta_study_enabled) {
        log.info "Study-level meta-analysis disabled (set --meta_study \"combined,cases,controls,meta\" to enable)"
    }
    if (meta_study_enabled) {
        def meta_study_list = meta_study_raw.toString()
            .split(',')
            .collect { it.trim().toLowerCase() }
            .findAll { it }
        def meta_study_set = meta_study_list as Set

        ch_meta_candidates = POST_GWAS_QC.out
            .map { sid, grp, prot, gwas, metrics, fig -> tuple(sid, grp, prot, gwas) }
            .mix(ch_meta_group_for_meta)
            .filter { sid, grp, prot, gwas -> meta_study_set.contains(grp.toString().toLowerCase()) }

        // Group all studies together for a given (group, protein)
        ch_for_meta = ch_meta_candidates
            .map { sid, grp, prot, gwas -> tuple("${grp}_${prot}", grp, prot, gwas) }
            .groupTuple(by: 0)
            .map { key, grps, prots, ss_list -> tuple(grps[0], prots[0], ss_list) }
            .combine(ch_studies_count)
            .filter { grp, prot, ss_list, count -> count > 1 && ss_list.size() == count }
            .map { grp, prot, ss_list, count -> tuple(grp, prot, ss_list) }

        CROSS_STUDY_META(ch_for_meta)
    }

    // 10. Post-QC on study-level meta outputs (METAL -> QC)
    def ch_meta_study_qc_out = Channel.empty()
    if (meta_study_enabled) {
        def ch_cross_meta_qc = CROSS_STUDY_META.out.map { grp, prot, meta_tbl, meta_info, meta_log ->
            tuple("meta-study", grp, prot, meta_tbl)
        }
        META_QC_STUDY(ch_cross_meta_qc)
        ch_meta_study_qc_out = META_QC_STUDY.out
    }
    
    // 14. Cohort Summary
    def ch_metrics_all = POST_GWAS_QC.out.map{sid, grp, prot, gwas, metrics, fig -> metrics.toString()}
    if (meta_group_enabled) {
        ch_metrics_all = ch_metrics_all.mix(ch_meta_group_qc_out.map { sid, grp, prot, gwas, metrics, fig -> metrics.toString() })
    }
    if (meta_study_enabled) {
        ch_metrics_all = ch_metrics_all.mix(ch_meta_study_qc_out.map { sid, grp, prot, gwas, metrics, fig -> metrics.toString() })
    }
    COHORT_SUMMARY(ch_metrics_all.collect())
}

#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("--phenotype_file"), type="character"),
  make_option(c("--sample_file"), type="character"),
  make_option(c("--study_id"), type="character"),
  make_option(c("--group"), type="character"),
  make_option(c("--group_column"), type="character", default=""),
  make_option(c("--cases_value"), type="character", default=""),
  make_option(c("--covariate_file"), type="character"),
  make_option(c("--covariates"), type="character", default=""),
  make_option(c("--include_proteins"), type="character", default=""),
  make_option(c("--outdir"), type="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(opt$outdir, "prepare_phenotypes.log")

# Helper to read protein list
get_protein_list <- function(input) {
  if (input == "") return(NULL)
  if (file.exists(input)) {
    prots <- readLines(input, warn = FALSE)
    prots <- trimws(prots)
    prots <- prots[prots != ""]
    return(prots)
  }
  return(trimws(unlist(strsplit(input, ","))))
}

standardize_ids <- function(df) {
  cols <- colnames(df)
  target <- grep("^#IID$|^IID$|^sampleid$|^sample_id$|^SampleID$", cols, ignore.case=TRUE, value=TRUE)
  if (length(target) > 0) {
    data.table::setnames(df, target[1], "IID")
  }
  return(df)
}

# 1. Load data
pheno <- data.table::fread(opt$phenotype_file)
pheno <- standardize_ids(pheno)

if (grepl("\\.psam$", opt$sample_file)) {
    samples <- data.table::fread(opt$sample_file, skip = "#IID")
} else if (grepl("\\.fam$|\\.sam$", opt$sample_file)) {
    samples <- data.table::fread(opt$sample_file, header = FALSE)
    data.table::setnames(samples, c("FID", "IID", "PID", "MID", "SEX", "PHENO"))
} else {
    samples <- data.table::fread(opt$sample_file)
}
samples <- standardize_ids(samples)

pheno[, IID := as.character(IID)]
samples[, IID := as.character(IID)]

common_samples <- intersect(pheno$IID, samples$IID)
if (length(common_samples) == 0) stop("Zero overlap found.")

pheno <- pheno[IID %in% common_samples]

# MATCH GENETIC DATA FORMAT: Set FID to 0
pheno[, FID := "0"]
data.table::setcolorder(pheno, c("FID", "IID"))

# 2. Determine columns to keep
all_cols <- colnames(pheno)
id_cols <- c("FID", "IID")
exclude_cols <- c("STUDY", "FID", "IID", "FID_1", "IID_1") # Common metadata

# Also exclude covariates provided by user
if (opt$covariates != "") {
  user_covs <- trimws(unlist(strsplit(opt$covariates, ",")))
  exclude_cols <- unique(c(exclude_cols, user_covs))
}

req_prots <- get_protein_list(opt$include_proteins)

if (!is.null(req_prots)) {
  # If proteins are provided, keep IDs + requested proteins
  # (Still excluding anything that might be a covariate by accident)
  valid_prots <- setdiff(intersect(req_prots, all_cols), exclude_cols)
  pheno <- pheno[, c(id_cols, valid_prots), with=FALSE]
} else {
  # If NO proteins provided, keep IDs + all numeric columns (excluding metadata/covariates)
  numeric_cols <- names(which(sapply(pheno, is.numeric)))
  trait_cols <- setdiff(numeric_cols, exclude_cols)
  pheno <- pheno[, c(id_cols, trait_cols), with=FALSE]
}

# 3. Group filtering
if (!is.null(opt$group) && opt$group != "combined") {
  covs <- data.table::fread(opt$covariate_file)
  covs <- standardize_ids(covs)
  covs[, IID := as.character(IID)]

  group_df <- NULL
  if (!is.null(opt$group_column) && opt$group_column != "" && opt$group_column %in% colnames(samples)) {
    group_df <- samples[, c("IID", opt$group_column), with=FALSE]
  } else if (!is.null(opt$group_column) && opt$group_column != "" && opt$group_column %in% colnames(covs)) {
    group_df <- covs[, c("IID", opt$group_column), with=FALSE]
  } else {
    stop(paste("Group column not found in sample or covariate file:", opt$group_column))
  }

  pheno <- merge(pheno, group_df, by="IID")
  if (opt$group == "cases") pheno <- pheno[get(opt$group_column) == opt$cases_value]
  else if (opt$group == "controls") pheno <- pheno[get(opt$group_column) != opt$cases_value]
  pheno[, (opt$group_column) := NULL]
  data.table::setcolorder(pheno, c("FID", "IID"))
}

# 4. Finalise
# Using NA for quantitative traits is standard for REGENIE
# No manual replacement with -9

full_out <- file.path(opt$outdir, "full.pheno")
data.table::fwrite(pheno, full_out, sep="\t", quote=FALSE, na="NA")

keep_out <- file.path(opt$outdir, "keep_samples.txt")
data.table::fwrite(pheno[, .(FID, IID)], keep_out, sep="\t", col.names=FALSE, quote=FALSE)

# 5. Summary statistics for proteins
prot_cols <- setdiff(colnames(pheno), c("FID", "IID"))
summary_list <- list()
for (p in prot_cols) {
  vals <- pheno[[p]]
  summary_list[[p]] <- data.frame(
    protein = p,
    n_total = length(vals),
    n_missing = sum(is.na(vals)),
    mean = mean(vals, na.rm=TRUE),
    sd = sd(vals, na.rm=TRUE),
    min = min(vals, na.rm=TRUE),
    max = max(vals, na.rm=TRUE)
  )
}
summary_df <- do.call(rbind, summary_list)
data.table::fwrite(summary_df, file.path(opt$outdir, "protein_summary.tsv"), sep="\t", quote=FALSE)

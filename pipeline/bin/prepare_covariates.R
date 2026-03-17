#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("--covariate_file"), type="character"),
  make_option(c("--sample_file"), type="character"),
  make_option(c("--study_id"), type="character"),
  make_option(c("--group"), type="character"),
  make_option(c("--group_column"), type="character", default=""),
  make_option(c("--cases_value"), type="character", default=""),
  make_option(c("--include_covariates"), type="character", default=""),
  make_option(c("--include_proteins"), type="character", default=""),
  make_option(c("--outdir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

standardize_ids <- function(df) {
  cols <- colnames(df)
  target <- grep("^#IID$|^IID$|^sampleid$|^sample_id$|^SampleID$", cols, ignore.case=TRUE, value=TRUE)
  if (length(target) > 0) {
    data.table::setnames(df, target[1], "IID")
  }
  return(df)
}

# 1. Load data
covs <- data.table::fread(opt$covariate_file)
covs <- standardize_ids(covs)

if (grepl("\\.psam$", opt$sample_file)) {
    samples <- data.table::fread(opt$sample_file, skip = "#IID")
} else if (grepl("\\.fam$|\\.sam$", opt$sample_file)) {
    samples <- data.table::fread(opt$sample_file, header = FALSE)
    data.table::setnames(samples, c("FID", "IID", "PID", "MID", "SEX", "PHENO"))
} else {
    samples <- data.table::fread(opt$sample_file)
}
samples <- standardize_ids(samples)

covs[, IID := as.character(IID)]
samples[, IID := as.character(IID)]

common_samples <- intersect(covs$IID, samples$IID)
if (length(common_samples) == 0) stop("Zero overlap found.")

covs <- covs[IID %in% common_samples]
# MATCH GENETIC DATA FORMAT: Set FID to 0
covs[, FID := "0"]
data.table::setcolorder(covs, c("FID", "IID"))

# Helper to read protein list (shared logic)
get_protein_list <- function(input) {
  if (input == "" || is.null(input)) return(NULL)
  if (file.exists(input)) {
    prots <- readLines(input, warn = FALSE)
    prots <- trimws(prots)
    prots <- prots[prots != ""]
    return(prots)
  }
  return(trimws(unlist(strsplit(input, ","))))
}

# 2. Covariate selection
all_cols <- colnames(covs)
id_cols <- c("FID", "IID")

if (opt$include_covariates != "") {
    req_covs <- trimws(unlist(strsplit(opt$include_covariates, ",")))
    if (!is.null(opt$group) && opt$group == "combined" && opt$group_column != "" &&
        !(opt$group_column %in% req_covs)) {
      req_covs <- c(req_covs, opt$group_column)
    }
    cols_to_keep <- c(id_cols, req_covs)
    covs <- covs[, intersect(cols_to_keep, all_cols), with=FALSE]
} else {
    # If NO explicit covariates provided, use all EXCEPT IDs and PROTEINS
    exclude_cols <- c("STUDY", "FID", "IID", "FID_1", "IID_1", "PlateId")
    
    # Also exclude proteins if a list or file is provided
    prots_to_exclude <- get_protein_list(opt$include_proteins)
    if (!is.null(prots_to_exclude)) {
        exclude_cols <- unique(c(exclude_cols, prots_to_exclude))
    }
    
    # For safety, only keep numeric columns as covariates if we are auto-selecting
    numeric_cols <- names(which(sapply(covs, is.numeric)))
    keep_cols <- setdiff(numeric_cols, exclude_cols)
    covs <- covs[, c(id_cols, keep_cols), with=FALSE]
}

# 3. Group filter
if (!is.null(opt$group) && opt$group != "combined") {
  group_df <- NULL
  if (!is.null(opt$group_column) && opt$group_column != "" && opt$group_column %in% colnames(samples)) {
    group_df <- samples[, c("IID", opt$group_column), with=FALSE]
  } else if (!is.null(opt$group_column) && opt$group_column != "" && opt$group_column %in% colnames(covs)) {
    group_df <- covs[, c("IID", opt$group_column), with=FALSE]
  } else {
    stop(paste("Group column not found in sample or covariate file:", opt$group_column))
  }

  covs <- merge(covs, group_df, by="IID")
  if (opt$group == "cases") covs <- covs[get(opt$group_column) == opt$cases_value]
  else if (opt$group == "controls") covs <- covs[get(opt$group_column) != opt$cases_value]

  if (!is.null(opt$group_column) && opt$group_column %in% colnames(covs)) {
    covs[, (opt$group_column) := NULL]
  }
}

# For combined group, include PHENO (group_column) as covariate
if (!is.null(opt$group) && opt$group == "combined" && !is.null(opt$group_column) && opt$group_column != "") {
  if (opt$group_column %in% colnames(samples)) {
    group_df <- samples[, c("IID", opt$group_column), with=FALSE]
    covs <- merge(covs, group_df, by="IID", all.x=TRUE)
  } else if (opt$group_column %in% colnames(covs)) {
    # already present
  } else {
    stop(paste("Group column not found in sample or covariate file:", opt$group_column))
  }
}

# 4. Process Columns
cov_cols <- setdiff(colnames(covs), c("FID", "IID"))
for (col in cov_cols) {
  if (all(is.na(covs[[col]]))) {
    cat("WARNING: Dropping all-NA covariate:", col, "\n")
    covs[, (col) := NULL]; next
  }

  unique_vals <- unique(na.omit(covs[[col]]))
  unique_count <- length(unique_vals)
  
  if (unique_count < 2) {
    cat("WARNING: Dropping invariant covariate:", col, "\n")
    covs[, (col) := NULL]; next
  }

  if (is.numeric(covs[[col]])) {
    # If it's a binary variable (e.g. 1/2), convert to 0/1 for REGENIE
    if (unique_count == 2) {
       sorted_vals <- sort(unique_vals)
       cat(">>> Info: Converting binary covariate", col, "(", sorted_vals[1], ",", sorted_vals[2], ") to 0/1\n")
       covs[, (col) := ifelse(get(col) == sorted_vals[1], 0, 1)]
    } else {
       # Continuous: Impute missing with mean but DON'T scale
       val_mean <- mean(covs[[col]], na.rm=TRUE)
       covs[is.na(get(col)), (col) := val_mean]
    }
  } else {
    # Categorical
    if (unique_count > (nrow(covs) * 0.5)) {
       cat("WARNING: Dropping covariate with too many levels:", col, "\n")
       covs[, (col) := NULL]; next
    }
    mode_val <- names(sort(table(covs[[col]]), decreasing = TRUE))[1]
    covs[is.na(get(col)), (col) := mode_val]
    
    mm <- model.matrix(as.formula(paste0("~ ", col)), data = covs)
    mm <- mm[, -1, drop=FALSE]
    covs <- cbind(covs, as.data.frame(mm))
    covs[, (col) := NULL]
  }
}

data.table::setcolorder(covs, c("FID", "IID"))
cov_file <- file.path(opt$outdir, "covariates.cov")
data.table::fwrite(covs, cov_file, sep="\t", quote=FALSE, na="NA")

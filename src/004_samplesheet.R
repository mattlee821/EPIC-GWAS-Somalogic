#!/usr/bin/env Rscript
# Script: src/004_samplesheet.R
# Purpose: Automatically generate the samplesheet.csv for the Nextflow GWAS pipeline
# based on the layout of genetics datasets and the main phenotype file.

# Adhere to USER coding rules: explicit package calls, use base R |>

library(optparse)

load_dot_env <- function(env_file = ".env") {
  if (file.exists(env_file)) {
    lines <- readLines(env_file, warn = FALSE)
    lines <- lines[!grepl("^\\s*(#|$)", lines)]
    for (line in lines) {
      parts <- strsplit(line, "=")[[1]]
      if (length(parts) >= 2) {
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = "="))
        value <- gsub("^['\"]|['\"]$", "", value)
        env_list <- list(value)
        names(env_list) <- key
        do.call(Sys.setenv, env_list)
      }
    }
  }
}

project_root <- Sys.getenv("GWAS_PROJECT_ROOT")
search_paths <- c(
  if (project_root != "") project_root else NULL,
  ".",
  "..",
  "../.."
)

for (p in search_paths) {
  env_p <- file.path(p, ".env")
  if (file.exists(env_p)) {
    load_dot_env(env_p)
    break
  }
}

get_env <- function(key, default = NULL) {
  val <- Sys.getenv(key)
  if (val == "") return(default)
  return(val)
}

read_sample_header <- function(sample_file) {
  hdr_lines <- readLines(sample_file, warn = FALSE)
  hdr_lines <- hdr_lines[!grepl("^\\s*$", hdr_lines)]
  hdr_lines <- hdr_lines[!grepl("^##", hdr_lines)]

  if (length(hdr_lines) == 0) {
    return(character())
  }

  hdr <- gsub("^#", "", hdr_lines[1])
  strsplit(trimws(hdr), "\\s+")[[1]]
}

option_list <- list(
  make_option(c("--pheno_file"), type="character", 
              default=get_env("GWAS_PHENO_FILE", "data/phenotype/phenofile.txt"),
              help="Path to the master phenotype text file"),
  make_option(c("--genetics_dir"), type="character", 
              default=get_env("GWAS_GENETICS_DIR", "data/genetics"),
              help="Path to the directory containing study genetic folders"),
  make_option(c("--out_file"), type="character", 
              default=get_env("GWAS_SAMPLESHEET", "pipeline/samplesheet.csv"),
              help="Output path for the samplesheet.csv"),
  make_option(c("--group_col"), type="character", default="PHENO",
              help="Column name for split analyses (e.g., case/control status). Default uses PHENO from .fam/.psam."),
  make_option(c("--cases_val"), type="character", default="2",
              help="Value in group_col mapping to 'cases' (default: 2)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (!file.exists(opt$pheno_file)) {
    stop("Phenotype file not found: ", opt$pheno_file)
}

if (!dir.exists(opt$genetics_dir)) {
    stop("Genetics directory not found: ", opt$genetics_dir)
}

# 1. Load the phenotype file to identify studies present
pheno <- data.table::fread(opt$pheno_file)

if (!"STUDY" %in% names(pheno)) {
    stop("Column 'STUDY' not found in phenotype file.")
}

unique_studies <- unique(pheno$STUDY)
cat("Found studies in phenotype file:\n")
print(unique_studies)

# 2. Iterate through studies and check for available genetic data
samplesheet_rows <- list()

for (study in unique_studies) {
    study_dir <- file.path(opt$genetics_dir, study)
    
    if (dir.exists(study_dir)) {
        # Search for pgen or bed files within the study
        pgen_files <- list.files(study_dir, pattern="\\.pgen$", full.names=TRUE, recursive=TRUE)
        bed_files <- list.files(study_dir, pattern="\\.bed$", full.names=TRUE, recursive=TRUE)
        
        genetic_hits <- if (length(pgen_files) > 0) pgen_files else bed_files
        
        if (length(genetic_hits) > 0) {
            # Use the first hit
            hit_path <- genetic_hits[1]
            ext_pattern <- if (length(pgen_files) > 0) "\\.pgen$" else "\\.bed$"
            pfile_prefix <- sub(ext_pattern, "", hit_path)
            
            # Check for BGEN directory
            bgen_dir <- file.path(study_dir, "bgen")
            if (!dir.exists(bgen_dir)) {
                bgen_dir <- file.path(dirname(pfile_prefix), "bgen")
            }
            
            # The sample file is either .psam or .fam/.sam
            sample_file <- if (length(pgen_files) > 0) paste0(pfile_prefix, ".psam") else paste0(pfile_prefix, ".fam")
            if (!file.exists(sample_file)) {
                fallback_fam <- paste0(pfile_prefix, ".fam")
                if (file.exists(fallback_fam)) {
                    sample_file <- fallback_fam
                }
            }
            
            # Function to ensure absolute paths even for non-existent files/prefixes
            get_abs_path <- function(p) {
                # If already absolute, return as is
                if (grepl("^/", p)) return(p)
                # Normalize against current working directory
                return(normalizePath(file.path(getwd(), p), mustWork = FALSE))
            }

            # Determine group column for case/control split
            group_col_use <- opt$group_col
            if (group_col_use != "") {
                if (grepl("\\.fam$|\\.sam$", sample_file)) {
                    # FAM/SAM has PHENO as 6th column by definition
                    if (group_col_use != "PHENO") {
                        cat("WARNING: sample_file is .fam/.sam but group_col is not PHENO. Using PHENO.\n")
                        group_col_use <- "PHENO"
                    }
                } else if (grepl("\\.psam$", sample_file)) {
                    if (file.exists(sample_file)) {
                        hdr_cols <- read_sample_header(sample_file)
                        group_col_found <- group_col_use %in% hdr_cols

                        if (!group_col_found && identical(group_col_use, "PHENO") && "PHENO1" %in% hdr_cols) {
                            group_col_found <- TRUE
                        }

                        if (!group_col_found) {
                            cat("WARNING: group_col '", group_col_use, "' not found in ", sample_file, ". Disabling split.\n", sep="")
                            group_col_use <- ""
                        }
                    } else {
                        cat("WARNING: sample_file not found to validate group_col. Disabling split.\n")
                        group_col_use <- ""
                    }
                }
            }

            # Record entry
            samplesheet_rows[[study]] <- data.frame(
              study_id = study,
              plink_bfile = get_abs_path(pfile_prefix),
              bgen_dir = get_abs_path(bgen_dir),
              sample_file = get_abs_path(sample_file),
              group_column = group_col_use,
              cases_value = opt$cases_val,
              stringsAsFactors = FALSE
            )
            cat("Matched genetic data for study:", study, " (Format:", if(length(pgen_files)>0) "PGEN" else "BED", ")\n")
        } else {
            cat("WARNING: No .pgen or .bed file found for study", study, "in", study_dir, "\n")
        }
    } else {
        cat("WARNING: Directory not found for study:", study_dir, "\n")
    }
}

# 3. Create flat data frame and save
if (length(samplesheet_rows) == 0) {
    stop("No matching genetic datasets found for any study. Cannot generate samplesheet.")
}

final_df <- do.call(rbind, samplesheet_rows)
rownames(final_df) <- NULL

# Ensure the output directory exists
out_dir <- dirname(opt$out_file)
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

data.table::fwrite(final_df, file=opt$out_file, sep=",", quote=FALSE, row.names=FALSE)
cat("\nSuccessfully wrote samplesheet to:", opt$out_file, "\n")

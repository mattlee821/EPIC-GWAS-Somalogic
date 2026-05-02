# Script: src/003_data-phenotype.R
# Purpose: Prepare phenotype inputs for the GWAS pipeline.

# environment ====
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

subset_seed <- as.integer(get_env("GWAS_PHENO_SUBSET_SEED", "20260421"))
set.seed(subset_seed)

# data somalogic ====
data_soma_path <- get_env("GWAS_EXTERNAL_SOMA_RDS")
if (is.null(data_soma_path) || !file.exists(data_soma_path)) {
  stop("GWAS_EXTERNAL_SOMA_RDS not found. Check your .env file.")
}
data_soma <- readRDS(data_soma_path)
data_features <- data_soma$data
meta_sample <- data_soma$others
meta_feature <- data_soma$aux

# format 
meta_sample <- meta_sample |>
  dplyr::mutate(IID = Idepic,
                IID_clean = gsub("_+", "_", Idepic))
data_features <- data_features |>
  dplyr::left_join(meta_sample |>
                     dplyr::select(SubjectID, IID, IID_clean),
                   by = "SubjectID") |>
  dplyr::select(SubjectID, IID, IID_clean, dplyr::everything()) |>
  dplyr::rename_with(
    ~ gsub("\\.", "-", gsub("^seq\\.", "", .x)),
    dplyr::starts_with("seq."))

protein_cols <- unique(meta_feature$SeqId)

# data genetic ====
depot_genetics_path <- get_env("GWAS_DEPOT_GENETICS_DIR")
if (is.null(depot_genetics_path) || !dir.exists(depot_genetics_path)) {
  stop("GWAS_DEPOT_GENETICS_DIR not found. Check your .env file.")
}

data_genetic <- list.files(path = depot_genetics_path, 
                           pattern = ".psam", all.files = T, full.names = T, recursive = T) |>
  purrr::map_dfr(function(f) {
    STUDY_name <- basename(f) |>
      tools::file_path_sans_ext()
    data.table::fread(f) |>
      dplyr::rename(IID = '#IID') |>
      dplyr::mutate(IID_clean = gsub("_+", "_", IID)) |>
      dplyr::mutate(STUDY = STUDY_name) |>
      dplyr::select(-SEX)}) |>
  dplyr::filter(IID %in% meta_sample$IID)

genetics_dir <- get_env("GWAS_GENETICS_DIR", "data/genetics")
data_genetic_neuro <- data.table::fread(file.path(genetics_dir, "Neuro_01/EPIC_GSA_imputed_QC_plink.psam")) |>
  dplyr::rename(IID = "#IID") |>
  dplyr::mutate(IID_clean = IID)  |>
  dplyr::select(-PHENO, -SEX) |>
  dplyr::mutate(STUDY = "Neuro_01") |>
  dplyr::filter(IID_clean %in% meta_sample$IID_clean)

data_genetic <- dplyr::bind_rows(
  data_genetic_neuro,
  data_genetic) |>
  dplyr::distinct(IID_clean, .keep_all = TRUE)

# make phenofile ====
data_pheno <- data_genetic |>
  dplyr::left_join(meta_sample |>
                     dplyr::select(IID_clean, Sex, Age_Recr, PlateId),
                   by = "IID_clean") |>
  dplyr::left_join(data_features |>
                     dplyr::select(IID_clean, dplyr::all_of(protein_cols)),
                   by = "IID_clean") |>
  dplyr::mutate(IID = IID_clean) |>
  dplyr::select(-IID_clean)

# write ====
if (!dir.exists("data/phenotype")) {dir.create("data/phenotype", recursive = TRUE)}

random_proteins <- sample(protein_cols, length(protein_cols))

write_trait_list <- function(filename, n_traits) {
  trait_subset <- random_proteins[seq_len(min(n_traits, length(random_proteins)))]
  writeLines(trait_subset, con = file.path("data/phenotype", filename))
  message("Wrote ", filename, " with ", length(trait_subset), " traits.")
}

data.table::fwrite(
  x = data_pheno,
  file = "data/phenotype/phenofile.txt",
  append = FALSE,
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
message("Wrote phenofile.txt with ", length(protein_cols), " proteins.")

write_trait_list("traits-10.txt", 10)
write_trait_list("traits-100.txt", 100)
write_trait_list("traits-1000.txt", 1000)

# environment ====
source("pipeline/bin/config.R")

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
                     dplyr::select(IID_clean, `10620-21`, `18878-15`, `4721-54`, `8323-163`),
                   by = "IID_clean") |>
  dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.))) |>
  dplyr::mutate(IID = IID_clean) |>
  dplyr::select(-IID_clean)

# write ====
if (!dir.exists("data/phenotype")) {dir.create("data/phenotype", recursive = TRUE)}
data.table::fwrite(x = data_pheno, file = "data/phenotype/phenofile-test.txt", 
                   append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

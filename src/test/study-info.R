# environment ====
source("pipeline/bin/config.R")

# data ====
data_soma_path <- get_env("GWAS_EXTERNAL_SOMA_RDS")
if (is.null(data_soma_path) || !file.exists(data_soma_path)) {
  stop("GWAS_EXTERNAL_SOMA_RDS not found. Check your .env file.")
}
data_soma <- readRDS(data_soma_path)
data_soma <- data_soma$others

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
      dplyr::mutate(STUDY = STUDY_name) |>
      dplyr::rename(IID = '#IID')}) |>
  dplyr::filter(IID %in% data_soma$Idepic) |>
  dplyr::select(-SEX) |>
  dplyr::left_join(data_soma |>
                     dplyr::select(Idepic, Sex),
                   by = c("IID" = "Idepic"))
data_genetic$IID <- gsub("____","_", data_genetic$IID)
data_genetic$IID <- gsub("__","_", data_genetic$IID)
data_genetic$IID <- gsub("__","_", data_genetic$IID)

epic4nd_archive <- get_env("GWAS_EPIC4ND_ARCHIVE")
if (is.null(epic4nd_archive) || !dir.exists(epic4nd_archive)) {
  stop("GWAS_EPIC4ND_ARCHIVE not found. Check your .env file.")
}
data_epic4nd <- data.table::fread(file.path(epic4nd_archive, "EPIC_GSA_imputed_QC_plink.psam")) |>
  dplyr::rename(IID = "#IID",
                Sex = SEX) |>
  dplyr::select(-PHENO) |>
  dplyr::mutate(STUDY = "EPIC4ND") 

# combine ====
# remove overlaps from data_genetic
epic4nd_unique <- data_epic4nd |>
  dplyr::anti_join(data_genetic, by = "IID")
data_genetic_f1 <- dplyr::bind_rows(
  data_genetic,
  epic4nd_unique
)
# remove overlaps from data_genetic
genetic_unique <- data_genetic |>
  dplyr::anti_join(data_epic4nd, by = "IID")
data_genetic_f2 <- dplyr::bind_rows(
  genetic_unique,
  data_epic4nd
)

# plot function ====
make_plot <- function(df, outfile) {
  
  data_plot <- df |>
    dplyr::mutate(Sex = haven::as_factor(Sex)) |>
    dplyr::count(STUDY, Sex) |>
    dplyr::group_by(STUDY) |>
    dplyr::mutate(total_n = sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      STUDY = forcats::fct_reorder(STUDY, total_n, .desc = TRUE)
    )
  
  totals_df <- data_plot |>
    dplyr::group_by(STUDY) |>
    dplyr::summarise(total_n = sum(n), .groups = "drop")
  
  N_100 <- totals_df |>
    dplyr::filter(total_n > 100) |>
    dplyr::summarise(overall_total = sum(total_n))
  
  first_small <- totals_df |>
    dplyr::filter(total_n < 100) |>
    dplyr::slice_min(order_by = STUDY) |>
    dplyr::pull(STUDY)
  
  vline_pos <- as.numeric(first_small) - 0.5
  
  p <- ggplot2::ggplot(data_plot, ggplot2::aes(x = STUDY, y = n, fill = Sex)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(
      xintercept = vline_pos,
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = n),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3
    ) +
    ggplot2::geom_text(
      data = totals_df,
      ggplot2::aes(x = STUDY, y = total_n, label = total_n),
      inherit.aes = FALSE,
      vjust = -0.4,
      size = 3.5
    ) +
    ggplot2::labs(
      subtitle = paste0(
        "Total N = ", sum(totals_df$total_n),
        "\nN > 100 = ", N_100$overall_total
      ),
      x = "STUDY",
      y = "Count",
      fill = "Sex"
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::expand_limits(y = max(totals_df$total_n) * 1.05)
  
  ggplot2::ggsave(
    outfile,
    plot = p,
    device = "tiff",
    bg = "white",
    width = 12,
    height = 10,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
  return(p)
}

# plot ====
analysis_dir <- get_env("GWAS_ANALYSIS_DIR", "analysis")
out_subdir <- file.path(analysis_dir, "001_data")
if (!dir.exists(out_subdir)) dir.create(out_subdir, recursive = TRUE)

p <- make_plot(
  data_genetic_f1,
  file.path(out_subdir, "study-info.tiff"))

p1 <- make_plot(
  data_genetic_f2,
  file.path(out_subdir, "study-info_alt.tiff"))

p2 <- cowplot::plot_grid(
  p + theme(legend.position = "none"),
  p1 + theme(legend.position = "bottom"),
  ncol = 1, labels = c("removing overlap from EPIC4ND",
                       "removing overlap from studies"), label_x = 0.2)
p2
ggplot2::ggsave(
  file.path(out_subdir, "study-info_combined.tiff"),
  plot = p2,
  device = "tiff",
  bg = "white",
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

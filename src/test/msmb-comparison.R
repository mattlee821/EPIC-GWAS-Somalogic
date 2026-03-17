rm(list = ls())
set.seed(821)

# environment ====
library(data.table)
library(ggplot2)
library(ggrepel)
source("pipeline/bin/config.R")

# source ====
BASE_DIR <- get_env("GWAS_ANALYSIS_DIR", "analysis")
STUDY_ID <- "Neuro_01"
PROT_ID <- "10620-21" # MSMB
GROUPS <- c("combined", "cases", "controls", "meta")
CHR_VAL <- 10

# data ====
load_gwas <- function(group) {
  # Standard results use 009_QC, within-study meta uses QC
  subdir <- if (group == "meta") "QC" else "009_QC"
  path <- file.path(BASE_DIR, STUDY_ID, PROT_ID, group, subdir, "gwas.tsv.gz")

  if (!file.exists(path) && group == "meta") {
    # Fallback to cross-study meta path if within-study meta isn't found
    path_cross <- file.path(BASE_DIR, "meta", PROT_ID, "meta", "QC", "gwas.tsv.gz")
    if (file.exists(path_cross)) path <- path_cross
  }

  if (!file.exists(path)) {
    message("Warning: File not found: ", path)
    return(NULL)
  }

  dt <- data.table::fread(path)
  # Subset to Chr 10 immediately
  dt <- dt[get("CHR") == CHR_VAL]
  dt[, "analysis" := group]

  return(dt)
}

# combine ====
datalist <- list()
for (grp in GROUPS) {
  datalist[[grp]] <- load_gwas(grp)
}
all_data <- data.table::rbindlist(datalist)

if (nrow(all_data) == 0) {
  stop("No data loaded. Check paths and group names.")
}
table(all_data$analysis)

# top 10 table ====
df <- all_data[, .SD[order(get("P"))][1:10], by = "analysis"] |>
  dplyr::arrange(ID)

# Save the table
data.table::fwrite(df, "temp/msmb-chr10.txt", sep = "\t")

# --- 3. Manhattan Plot with Facet Wrap ---
# Prepare plot data
group_order <- c("combined", "meta", "controls", "cases")

min <- df |>
  dplyr::slice_min(order_by = P, n = 1, with_ties = FALSE) |>
  dplyr::pull(POS)
up <- min + 100000
down <- min - 100000

df_plot <- all_data |>
  dplyr::mutate(
    log10P = -log10(P),
    analysis = factor(analysis, levels = group_order)
  ) |>
  dplyr::filter(POS < up & POS > down)

# Prepare labels (deduplicate IDs to avoid over-labeling same SNP across groups if desired,
# but user asked to label points across each study/analysis)
label_dt <- df |>
  dplyr::mutate(
    log10P = -log10(P),
    analysis = factor(analysis, levels = group_order)
  )

# Create the plot
p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = POS, y = log10P)) +
  # Background points
  ggplot2::geom_point(alpha = 0.4, color = "grey50", size = 0.8) +
  # Highlight top points
  ggplot2::geom_point(data = label_dt, ggplot2::aes(color = analysis), size = 1.2) +
  # Labels
  ggrepel::geom_label_repel(
    data = label_dt,
    ggplot2::aes(label = ID),
    size = 2.5, force = 10,
    max.overlaps = 20, min.segment.length = 0,
    box.padding = 0.5
  ) +
  # Faceting - one row per analysis
  ggplot2::facet_wrap(~analysis, ncol = 1, scales = "free_y") +
  # Style
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = paste("MSMB (10620-21)"),
    x = paste("CHR", CHR_VAL, "(bp)"),
    y = "-log10(P-value)",
    color = ""
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey90"),
    strip.text = ggplot2::element_text(face = "bold")
  ) +
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "none") +
  # Significance lines
  ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "solid", alpha = 0.5)

# Save the plot
ggplot2::ggsave("temp/msmb-chr10.png", p, width = 12, height = 16, dpi = 300)

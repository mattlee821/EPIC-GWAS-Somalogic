#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("--metrics_files"), type="character", help="Comma-separated list of metrics files"),
  make_option(c("--metrics_list"), type="character", help="File containing list of metrics files (one per line)"),
  make_option(c("--outdir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Parse paths from either comma-separated string or list file
metrics_paths <- character(0)

if (!is.null(opt$metrics_list) && file.exists(opt$metrics_list)) {
  metrics_paths <- readLines(opt$metrics_list) |> 
    trimws()
  metrics_paths <- metrics_paths[nchar(metrics_paths) > 0]
} else if (!is.null(opt$metrics_files)) {

  metrics_paths <- trimws(unlist(strsplit(opt$metrics_files, ",")))
}

# Combine all metrics
if (length(metrics_paths) > 0 && metrics_paths[1] != "") {
  metrics_list <- lapply(metrics_paths, data.table::fread)
  metrics_all <- data.table::rbindlist(metrics_list, fill=TRUE)

  drop_cols <- c("flag_raw", "flag_filtered", "info_score")
  present <- intersect(drop_cols, colnames(metrics_all))
  if (length(present) > 0) {
    # Using explicit data.table::set for column removal
    data.table::set(metrics_all, j = present, value = NULL)
  }

  metrics_out <- file.path(opt$outdir, "all_metrics.tsv")
  data.table::fwrite(metrics_all, metrics_out, sep="\t", quote=FALSE)

  # Lambda GC Distribution Plot
  lambda_col <- NULL
  if ("l_filt" %in% colnames(metrics_all)) {
    lambda_col <- "l_filt"
  } else if ("l_raw" %in% colnames(metrics_all)) {
    lambda_col <- "l_raw"
  } else if ("lambda_gc_filtered" %in% colnames(metrics_all)) {
    lambda_col <- "lambda_gc_filtered"
  } else if ("lambda_gc_raw" %in% colnames(metrics_all)) {
    lambda_col <- "lambda_gc_raw"
  }


  if (!is.null(lambda_col)) {
    library(ggplot2)
    p <- ggplot2::ggplot(metrics_all, ggplot2::aes(x=.data[[lambda_col]])) +
      ggplot2::geom_histogram(bins=50, fill="steelblue", color="black") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title="Distribution of Genomic Inflation (Lambda GC)",
           x=lambda_col, y="Count")
    ggplot2::ggsave(file.path(opt$outdir, "lambda_gc_distribution.png"), p, width=8, height=6)
  }
}


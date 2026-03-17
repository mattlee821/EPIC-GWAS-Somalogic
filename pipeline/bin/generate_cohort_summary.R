#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("--metrics_files"), type="character"),
  make_option(c("--outdir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Parse lists
metrics_paths <- trimws(unlist(strsplit(opt$metrics_files, ",")))

# Combine all metrics
if (length(metrics_paths) > 0 && metrics_paths[1] != "") {
  metrics_list <- lapply(metrics_paths, data.table::fread)
  metrics_all <- data.table::rbindlist(metrics_list, fill=TRUE)

  drop_cols <- c("flag_raw", "flag_filtered", "info_score")
  present <- intersect(drop_cols, colnames(metrics_all))
  if (length(present) > 0) {
    metrics_all[, (present) := NULL]
  }

  metrics_out <- file.path(opt$outdir, "all_metrics.tsv")
  data.table::fwrite(metrics_all, metrics_out, sep="\t", quote=FALSE)

  # Lambda GC Distribution Plot
  lambda_col <- NULL
  if ("lambda_gc_filtered" %in% colnames(metrics_all)) {
    lambda_col <- "lambda_gc_filtered"
  } else if ("lambda_gc_raw" %in% colnames(metrics_all)) {
    lambda_col <- "lambda_gc_raw"
  } else if ("lambda_gc" %in% colnames(metrics_all)) {
    lambda_col <- "lambda_gc"
  }

  if (!is.null(lambda_col)) {
    library(ggplot2)
    p <- ggplot(metrics_all, aes(x=.data[[lambda_col]])) +
      geom_histogram(bins=50, fill="steelblue", color="black") +
      theme_minimal() +
      labs(title="Distribution of Genomic Inflation (Lambda GC)",
           x=lambda_col, y="Count")
    ggsave(file.path(opt$outdir, "lambda-gc-distribution.png"), p, width=8, height=6)
  }
}

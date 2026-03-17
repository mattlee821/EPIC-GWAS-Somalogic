#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(ggplot2)

option_list = list(
  make_option(c("--cov_file"), type="character"),
  make_option(c("--pcs_file"), type="character"),
  make_option(c("--eval_file"), type="character"),
  make_option(c("--outdir"), type="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# 1. Load data
covs <- fread(opt$cov_file)
pcs <- fread(opt$pcs_file, header = FALSE)

# Format PC file
if (ncol(pcs) >= 12) {
  setnames(pcs, c("FID", "IID", paste0("PC", 1:(ncol(pcs)-2))))
} else {
  # If header exists, try reading with header
  pcs <- fread(opt$pcs_file)
  if (!all(c("FID", "IID") %in% colnames(pcs))) {
    stop("PC file missing FID/IID columns.")
  }
}

covs[, FID := as.character(FID)]
covs[, IID := as.character(IID)]
pcs[, FID := as.character(FID)]
pcs[, IID := as.character(IID)]

# 2. Merge
merged <- merge(covs, pcs, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)

out_file <- file.path(opt$outdir, "covariates.cov")
fwrite(merged, out_file, sep = "\t", quote = FALSE, na = "NA")

# 3. Scree Plot
if (file.exists(opt$eval_file)) {
  evals <- fread(opt$eval_file, header = FALSE)
  colnames(evals) <- "Val"
  evals[, PC := seq_len(.N)]
  evals[, PVE := (Val / sum(Val)) * 100]
  
  p <- ggplot(evals, aes(x = factor(PC), y = PVE)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", PVE)), vjust = -0.5, size = 3) +
    theme_minimal() +
    labs(title = "Scree Plot",
         x = "Principal Component",
         y = "Proportion of Variance Explained (%)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

  ggsave(file.path(opt$outdir, "scree_plot.png"), p, width = 8, height = 5, dpi = 300)
}

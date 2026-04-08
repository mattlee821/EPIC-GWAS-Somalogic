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
pcs <- fread(opt$pcs_file)

# Format PC file (plink2 pca.eigenvec often uses '#FID' as first header)
if ("#FID" %in% colnames(pcs) && !"FID" %in% colnames(pcs)) {
  setnames(pcs, "#FID", "FID")
}

# If headerless, fall back and create names
if (!all(c("FID", "IID") %in% colnames(pcs))) {
  pcs <- fread(opt$pcs_file, header = FALSE)
  if (ncol(pcs) < 3) {
    stop("PC file has too few columns (expected FID IID PC1 ...).")
  }
  setnames(pcs, c("FID", "IID", paste0("PC", 1:(ncol(pcs) - 2))))
}

covs[, FID := as.character(FID)]
covs[, IID := as.character(IID)]
pcs[, FID := as.character(FID)]
pcs[, IID := as.character(IID)]

# 2. Merge
merged <- merge(covs, pcs, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)

# Ensure covariate column names are safe for downstream tools (keep FID/IID intact)
id_cols <- c("FID", "IID")
other_cols <- setdiff(colnames(merged), id_cols)
if (length(other_cols) > 0) {
  safe_other <- make.names(other_cols, unique = TRUE)
  setnames(merged, old = other_cols, new = safe_other)
}

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

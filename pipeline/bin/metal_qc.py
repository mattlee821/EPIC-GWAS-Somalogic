#!/usr/bin/env python3
import argparse
import os
import numpy as np
import polars as pl
from scipy.stats import chi2

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True, help="Single METAL output file")
    parser.add_argument("--protein_id", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()

def calculate_lambda(p_values):
    p = p_values.drop_nulls()
    if len(p) == 0: return np.nan
    p_arr = np.clip(p.to_numpy(), 1e-300, 1.0)
    return np.nanmedian(chi2.isf(p_arr, df=1)) / chi2.ppf(0.5, df=1)

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mapping = {
        "MarkerName": "ID",
        "Allele1": "EA",
        "Allele2": "OA",
        "Freq1": "EAF",
        "Effect": "BETA",
        "StdErr": "SE",
        "P-value": "P",
        "Weight": "N",
        "Chromosome": "CHR",
        "Position": "POS",
        "CHROMOSOME": "CHR",
        "POSITION": "POS",
    }

    # Read the raw METAL output or an already-standardized table, then keep only
    # the columns needed for final output and QC metrics.
    df = pl.read_csv(args.input_file, separator="\t", ignore_errors=True)
    df = df.rename({old: new for old, new in mapping.items() if old in df.columns})
    keep_cols = ["CHR", "POS", "ID", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
    df = df.select([c for c in keep_cols if c in df.columns])

    cast_exprs = []
    for col in ("P", "N", "EAF", "BETA", "SE"):
        if col in df.columns:
            cast_exprs.append(pl.col(col).cast(pl.Float64, strict=False))
    if "POS" in df.columns:
        cast_exprs.append(pl.col("POS").cast(pl.Int64, strict=False))
    if cast_exprs:
        df = df.with_columns(cast_exprs)

    # Ensure CHR and POS columns are present and correctly handled
    needs_split = False
    if "CHR" not in df.columns or "POS" not in df.columns:
        needs_split = True
    elif df["CHR"].null_count() == len(df) or df["POS"].null_count() == len(df):
        needs_split = True

    if needs_split and "ID" in df.columns:
        id_clean = df["ID"].cast(pl.Utf8).str.replace(":", "_")
        id_parts = id_clean.str.split("_")

        new_cols = []
        if "CHR" not in df.columns:
            new_cols.append(id_parts.list.get(0).str.replace("(?i)chr", "").alias("CHR"))
        elif df["CHR"].null_count() > 0:
            # Fill nulls only
            new_cols.append(pl.when(pl.col("CHR").is_null())
                            .then(id_parts.list.get(0).str.replace("(?i)chr", ""))
                            .otherwise(pl.col("CHR").cast(pl.Utf8).str.replace("(?i)chr", ""))
                            .alias("CHR"))
        else:
             new_cols.append(pl.col("CHR").cast(pl.Utf8).str.replace("(?i)chr", "").alias("CHR"))

        if "POS" not in df.columns:
            new_cols.append(id_parts.list.get(1).cast(pl.Int64, strict=False).alias("POS"))
        elif df["POS"].null_count() > 0:
            new_cols.append(pl.when(pl.col("POS").is_null())
                            .then(id_parts.list.get(1).cast(pl.Int64, strict=False))
                            .otherwise(pl.col("POS"))
                            .alias("POS"))
        
        df = df.with_columns(new_cols)
    else:
        if "CHR" in df.columns:
            df = df.with_columns(pl.col("CHR").cast(pl.Utf8).str.replace("(?i)chr", ""))

    if "EA" in df.columns:
        df = df.with_columns(pl.col("EA").str.to_uppercase())
    if "OA" in df.columns:
        df = df.with_columns(pl.col("OA").str.to_uppercase())

    l_raw = calculate_lambda(df["P"])

    if "N" in df.columns:
        max_n = df["N"].max()
        if max_n is not None and max_n > 0:
            df_filt = df.filter(pl.col("N") >= 0.5 * max_n)
        else:
            df_filt = df
    else:
        df_filt = df
    
    l_filt = calculate_lambda(df_filt["P"])

    standard_cols = ["CHR", "POS", "ID", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]

    for col in standard_cols:
        if col not in df_filt.columns:
            df_filt = df_filt.with_columns(pl.lit(None).alias(col))

    final_df = df_filt.select(standard_cols)

    gwas_out = os.path.join(args.outdir, "meta.tsv.gz")
    final_df.lazy().sink_csv(gwas_out, separator="\t", compression="gzip")

    metrics_path = os.path.join(args.outdir, "metrics.tsv")
    metrics_df = pl.DataFrame([{
        "study": args.study,
        "group": args.group,
        "id": args.protein_id,
        "l_raw": l_raw,
        "l_filt": l_filt
    }])
    metrics_df.write_csv(metrics_path, separator="\t")

if __name__ == "__main__":
    main()

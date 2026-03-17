#!/usr/bin/env python3
import argparse
import os
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import chi2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import pyarrow as pa
    import pyarrow.csv as pacsv
    _HAS_PYARROW = True
except Exception:
    _HAS_PYARROW = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Finalize GWAS results (combine, QC, plot, export).")
    parser.add_argument("--input_files", required=True, help="Comma-separated list of REGENIE step2 files")
    parser.add_argument("--study", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--protein_id", required=True)
    parser.add_argument("--info_score", type=float, default=0.3)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


# Utility to get path from env if needed in the future
def get_path_env(key: str, default: str) -> str:
    return os.getenv(key, default)


def read_regenie_files(files: List[str]) -> pd.DataFrame:
    df = None
    for f in files:
        chunk = pd.read_csv(f, sep=r"\s+", engine="python", compression="infer")
        chunk.columns = [c.strip() for c in chunk.columns]
        if df is None:
            df = chunk
        else:
            df = pd.concat([df, chunk], ignore_index=True)
        del chunk
    if df is None:
        raise RuntimeError("No input files found for finalise.py")
    return df


def derive_chr_pos_from_id(df: pd.DataFrame, id_col: str = "ID") -> None:
    if id_col not in df.columns:
        return
    s = df[id_col].astype(str)
    # Try patterns like "15_33001734_..." or "15:33001734..."
    m = s.str.extract(r"^(?:chr)?(?P<chr>[0-9XYMT]+)[_:](?P<pos>[0-9]+)")
    if m["chr"].isna().all():
        m = s.str.extract(r"^(?:chr)?(?P<chr>[0-9XYMT]+)_(?P<pos>[0-9]+)")
    if "chr" in m.columns:
        df["CHR"] = m["chr"]
    if "pos" in m.columns:
        df["POS"] = pd.to_numeric(m["pos"], errors="coerce")


def add_chr_num(df: pd.DataFrame, chrom_col: str) -> None:
    df["CHR_NUM"] = pd.to_numeric(df[chrom_col], errors="coerce")
    mask_chr = df["CHR_NUM"].isna()
    if mask_chr.any():
        ch = df.loc[mask_chr, chrom_col].astype(str).str.replace("chr", "", case=False).str.upper()
        mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        df.loc[mask_chr, "CHR_NUM"] = ch.map(mapping)


def manhattan_plot(ax, df: pd.DataFrame, chr_col: str, pos_col: str, p_col: str, title: str) -> None:
    if df.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_axis_off()
        return

    chr_vals = df[chr_col]
    pos_vals = df[pos_col]
    pvals = df[p_col].clip(lower=1e-300)

    mask = chr_vals.notna() & pos_vals.notna() & pvals.notna()
    if not mask.any():
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_axis_off()
        return

    chr_max = pos_vals[mask].groupby(chr_vals[mask]).max().sort_index()
    offsets = chr_max.cumsum() - chr_max
    pos_cum = pos_vals + chr_vals.map(offsets)
    neglogp = -np.log10(pvals)

    colors = ["#4E79A7", "#F28E2B"]
    for i, chr_val in enumerate(chr_max.index.tolist()):
        chr_mask = mask & (chr_vals == chr_val)
        ax.scatter(pos_cum[chr_mask], neglogp[chr_mask], c=colors[i % 2], s=2, linewidths=0)

    midpoints = (offsets + chr_max / 2.0).to_dict()
    ax.set_xticks(list(midpoints.values()))
    ax.set_xticklabels([str(k) for k in midpoints.keys()], fontsize=7)
    ax.set_ylabel("-log10(P)")
    ax.set_title(title, fontsize=9)


def qq_plot(ax, pvals: pd.Series, title: str) -> None:
    pvals = pvals.dropna()
    pvals = pvals[pvals > 0]
    if pvals.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_axis_off()
        return

    n = pvals.shape[0]
    exp = -np.log10(np.linspace(1.0 / (n + 1), n / (n + 1), n))
    obs = -np.log10(np.sort(pvals.values))
    ax.scatter(exp, obs, s=2, color="#4E79A7")
    maxv = max(exp.max(), obs.max())
    ax.plot([0, maxv], [0, maxv], color="red", linewidth=1)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(title, fontsize=9)


def main() -> int:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    files = [f for f in (x.strip() for x in args.input_files.split(",")) if f]
    sumstats = read_regenie_files(files)

    # Harmonize common METAL-style column names
    if "ID" not in sumstats.columns and "MarkerName" in sumstats.columns:
        sumstats["ID"] = sumstats["MarkerName"]
    if "EA" not in sumstats.columns and "Allele1" in sumstats.columns:
        sumstats["EA"] = sumstats["Allele1"]
    if "OA" not in sumstats.columns and "Allele2" in sumstats.columns:
        sumstats["OA"] = sumstats["Allele2"]
    if "EAF" not in sumstats.columns and "Freq1" in sumstats.columns:
        sumstats["EAF"] = sumstats["Freq1"]
    if "BETA" not in sumstats.columns and "Effect" in sumstats.columns:
        sumstats["BETA"] = sumstats["Effect"]
    if "SE" not in sumstats.columns and "StdErr" in sumstats.columns:
        sumstats["SE"] = sumstats["StdErr"]
    if "P" not in sumstats.columns and "P-value" in sumstats.columns:
        sumstats["P"] = sumstats["P-value"]
    if "N" not in sumstats.columns and "Weight" in sumstats.columns:
        sumstats["N"] = sumstats["Weight"]

    # Strip quotes and normalize allele casing
    if "ID" in sumstats.columns:
        sumstats["ID"] = sumstats["ID"].astype(str).str.replace('"', "", regex=False).str.strip()
    for col in ["EA", "OA", "ALLELE0", "ALLELE1", "Allele1", "Allele2"]:
        if col in sumstats.columns:
            sumstats[col] = (
                sumstats[col]
                .astype(str)
                .str.replace('"', "", regex=False)
                .str.strip()
                .str.upper()
            )

    # Ensure we have a P-value column
    if "LOG10P" in sumstats.columns:
        sumstats["P"] = np.power(10.0, -sumstats["LOG10P"].astype(float))
    elif "P" in sumstats.columns:
        sumstats["P"] = pd.to_numeric(sumstats["P"], errors="coerce")
    else:
        # Try common P-value column variants
        for cand in ["PVAL", "PVALUE", "Pvalue", "P-value"]:
            if cand in sumstats.columns:
                sumstats["P"] = pd.to_numeric(sumstats[cand], errors="coerce")
                break
        if "P" not in sumstats.columns:
            raise RuntimeError(f"Missing P/LOG10P column in input sumstats. Columns: {list(sumstats.columns)}")

    # 1-4: load, format columns, compute P, keep only necessary columns
    keep_cols = [
        "CHROM",
        "GENPOS",
        "ID",
        "EA",
        "OA",
        "EAF",
        "ALLELE0",
        "ALLELE1",
        "A1FREQ",
        "BETA",
        "SE",
        "N",
        "INFO",
        "TEST",
        "CHISQ",
        "LOG10P",
        "P",
    ]
    keep_cols = [c for c in keep_cols if c in sumstats.columns]
    sumstats = sumstats[keep_cols]

    # Derive CHR/POS if not available (e.g., METAL outputs)
    if "CHROM" not in sumstats.columns and "CHR" not in sumstats.columns:
        derive_chr_pos_from_id(sumstats, "ID")

    chrom_col = "CHROM" if "CHROM" in sumstats.columns else ("CHR" if "CHR" in sumstats.columns else None)
    pos_col = "GENPOS" if "GENPOS" in sumstats.columns else ("POS" if "POS" in sumstats.columns else None)
    if chrom_col is not None:
        add_chr_num(sumstats, chrom_col)

    # 5: metrics on raw
    n_before = sumstats.shape[0]
    lambda_gc_raw = np.nan
    if n_before > 0:
        pvals_raw = sumstats["P"].clip(lower=1e-300)
        chisq_raw = chi2.isf(pvals_raw, df=1)
        lambda_gc_raw = np.nanmedian(chisq_raw) / chi2.ppf(0.5, df=1)

    # 6: raw figures
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=300)
    if pos_col is None:
        pos_col = "POS" if "POS" in sumstats.columns else "GENPOS"
    manhattan_plot(axes[0, 0], sumstats, "CHR_NUM", pos_col, "P", f"Pre-QC Manhattan: {args.protein_id}")
    qq_plot(axes[0, 1], sumstats["P"], "Pre-QC QQ")

    # 7: filter raw data
    mask = pd.Series(True, index=sumstats.index)
    if "INFO" in sumstats.columns:
        mask &= sumstats["INFO"] >= args.info_score
    if "N" in sumstats.columns:
        max_n = sumstats["N"].max(skipna=True)
        if pd.notna(max_n):
            mask &= sumstats["N"] >= 0.5 * max_n
    if "BETA" in sumstats.columns and "SE" in sumstats.columns:
        mask &= (sumstats["BETA"].abs() <= 10) & (sumstats["SE"] <= 10) & (sumstats["SE"] > 0)

    sumstats = sumstats.loc[mask].reset_index(drop=True)

    # 8: save filtered data (formatted output)
    rename_map = {
        "CHROM": "CHR",
        "GENPOS": "POS",
        "ALLELE0": "OA",
        "ALLELE1": "EA",
        "A1FREQ": "EAF",
    }
    sumstats.rename(columns=rename_map, inplace=True)

    if "CHISQ" not in sumstats.columns:
        sumstats["CHISQ"] = chi2.isf(sumstats["P"].clip(lower=1e-300), df=1)

    required_cols = ["CHR", "POS", "ID", "EA", "OA", "EAF", "BETA", "SE", "P", "N", "CHISQ"]
    for col in required_cols:
        if col not in sumstats.columns:
            sumstats[col] = np.nan

    # Keep CHR_NUM for plotting, drop later
    save_cols = [c for c in required_cols if c in sumstats.columns]
    
    out_gwas = os.path.join(args.outdir, "gwas.tsv.gz")
    if _HAS_PYARROW:
        table = pa.Table.from_pandas(sumstats[save_cols], preserve_index=False)
        write_opts = pacsv.WriteOptions(include_header=True, delimiter="\t")
        with pa.output_stream(out_gwas, compression="gzip") as sink:
            pacsv.write_csv(table, sink, write_options=write_opts)
    else:
        compression_opts = {"method": "gzip", "compresslevel": 1}
        sumstats.to_csv(out_gwas, sep="\t", index=False, columns=save_cols, compression=compression_opts)

    # 9: metrics on filtered data
    n_after = sumstats.shape[0]
    lambda_gc_filtered = np.nan
    if n_after > 0:
        pvals_f = sumstats["P"].clip(lower=1e-300)
        chisq_f = chi2.isf(pvals_f, df=1)
        lambda_gc_filtered = np.nanmedian(chisq_f) / chi2.ppf(0.5, df=1)

    # 10: save metrics
    metrics = pd.DataFrame([
        {
            "study": args.study,
            "group": args.group,
            "protein_id": args.protein_id,
            "n_variants_before": int(n_before),
            "n_variants_after": int(n_after),
            "lambda_gc_raw": float(lambda_gc_raw) if pd.notna(lambda_gc_raw) else np.nan,
            "lambda_gc_filtered": float(lambda_gc_filtered) if pd.notna(lambda_gc_filtered) else np.nan,
        }
    ])
    metrics_path = os.path.join(args.outdir, "metrics.tsv")
    metrics.to_csv(metrics_path, sep="\t", index=False)

    # 11-12: filtered figures + combined figure
    plot_chr_col = "CHR_NUM" if "CHR_NUM" in sumstats.columns else "CHR"
    manhattan_plot(axes[1, 0], sumstats, plot_chr_col, "POS", "P", f"Post-QC Manhattan: {args.protein_id}")
    qq_plot(axes[1, 1], sumstats["P"], "Post-QC QQ")
    fig.tight_layout()
    fig_path = os.path.join(args.outdir, "figure.png")
    fig.savefig(fig_path)
    plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

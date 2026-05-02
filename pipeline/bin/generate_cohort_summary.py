#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from config import read_table_auto


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics_files", default="")
    parser.add_argument("--metrics_list", default="")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    metrics_paths: list[str] = []
    if args.metrics_list and Path(args.metrics_list).exists():
        metrics_paths = [line.strip() for line in Path(args.metrics_list).read_text().splitlines() if line.strip()]
    elif args.metrics_files:
        metrics_paths = [part.strip() for part in args.metrics_files.split(",") if part.strip()]

    if not metrics_paths:
        return

    metrics_frames = [read_table_auto(path) for path in metrics_paths]
    metrics_all = pd.concat(metrics_frames, ignore_index=True, sort=False)

    drop_columns = ["flag_raw", "flag_filtered", "info_score"]
    present = [column for column in drop_columns if column in metrics_all.columns]
    if present:
        metrics_all = metrics_all.drop(columns=present)

    metrics_all.to_csv(outdir / "all_metrics.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()

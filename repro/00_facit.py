#!/usr/bin/env python3
"""Reads HPA's published per-gene log-rank p-values (the answer key) and tabulates, per cancer cohort, how many genes were tested,
how many are called prognostic (p < 0.001), and the median p. No new statistics: only HPA's own numbers.
Input: data/hpa/cancer_prognostic_data.tsv (HPA v25.1). Output: repro/ut/00_facit_per_cancer.csv"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "data" / "hpa" / "cancer_prognostic_data.tsv", sep="\t")
cols = list(d.columns[3:])
assert (d[cols].notna().sum(axis=1) <= 1).all(), "more than one label column filled"
d["p"] = d[cols].bfill(axis=1).iloc[:, 0]
g = d.dropna(subset=["p"]).groupby("Cancer").agg(genes_tested=("p", "size"), prognostic=("p", lambda x: int((x < 1e-3).sum())), median_p=("p", "median"))
g["share_prognostic"] = (g["prognostic"] / g["genes_tested"]).round(4)
g["expected_if_test_valid_and_null"] = (g["genes_tested"] * 1e-3).round(1)
g = g.sort_values("prognostic", ascending=False)
(ROOT / "repro" / "ut").mkdir(parents=True, exist_ok=True)
g.to_csv(ROOT / "repro" / "ut" / "00_facit_per_cancer.csv")
print(g.to_string())

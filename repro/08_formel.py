#!/usr/bin/env python3
"""The closed-form correction for a maximally selected log-rank statistic, laid next to the permutation calibration from 03.

Miller & Siegmund (Biometrics 1982), applied to the log-rank test by Lausen & Schumacher (Biometrics 1992). As implemented in the R
package maxstat (pLausen92; formula confirmed against that source code, see docs/04_LITTERATUR.md section 2; the 1992 paper itself is
behind a paywall and has not been read in the original):
    p_corr = phi(z) * (z - 1/z) * ln[ e2 (1 - e1) / ((1 - e2) e1) ] + 4 phi(z) / z,      z = Phi^-1(1 - p_min / 2)
with [e1, e2] the quantile range searched for the cutoff. HPA searches the 20th to 80th percentile, so the log term is ln 16.
The approximation is asymptotic and meant for small p_min; values above 1 are set to 1.

Reported: (1) what p_min = 0.001 really means according to the formula, and which p_min gives a real 0.001; (2) per cohort, the formula's
false-positive rate at HPA's threshold next to the permutation estimate from 03; (3) per gene, the ratio between the permutation-calibrated p and the
formula's p (both are monotone in p_min, so rank agreement is trivially perfect and is not reported), and how many genes pass Benjamini-Hochberg 5 % with each.
Input: repro/ut/03_<ABBR>_gener.csv, repro/ut/03_sammanfattning.csv. Output: repro/ut/08_formel.csv. No randomness."""
import importlib.util
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent; spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
UT = rep.UT


def p_lausen92(p_min, e1=0.20, e2=0.80):
    p_min = np.clip(np.asarray(p_min, float), 1e-300, 1); z = norm.isf(p_min / 2); z = np.maximum(z, 1e-12)
    p = norm.pdf(z) * (z - 1 / z) * np.log(e2 * (1 - e1) / ((1 - e2) * e1)) + 4 * norm.pdf(z) / z
    return np.where(z > 1, np.clip(p, 0, 1), 1.0)     # the approximation breaks down for z <= 1 (p_min >= 0.32)


def bh(p):
    n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]; ut = np.empty(n); ut[o] = np.minimum(q, 1); return ut


if __name__ == "__main__":
    a = float(p_lausen92(1e-3)); need = brentq(lambda x: float(p_lausen92(x)) - 1e-3, 1e-9, 1e-3)
    print(f"formula: p_min = 0.001 corresponds to p_corr = {a:.4f} ({a / 1e-3:.1f} x nominal); a real 0.001 needs p_min < {need:.2e}")
    s3 = pd.read_csv(UT / "03_sammanfattning.csv").set_index("cohort"); rows = []
    for ab in rep.NAMN:
        d = pd.read_csv(UT / f"03_{ab}_gener.csv").dropna(subset=["p_ours"]); d["p_formula"] = p_lausen92(d.p_ours.to_numpy()); d["q_formula"] = bh(d.p_formula.to_numpy()); pr = d.p_ours < 1e-3
        small = d[d.p_ours < 0.01]
        rows.append(dict(cohort=ab, events=int(s3.events[ab]), alpha_permutation=s3.alpha_eff[ab], alpha_formula=round(a, 5), ratio_perm_to_formula=round(s3.alpha_eff[ab] / a, 2),
                         median_ratio_pcal_to_formula_below_01=round(float(np.median(small.p_cal / small.p_formula)), 2) if len(small) else np.nan,
                         hpa_prognostic=int(pr.sum()), pass_bh05_permutation=int((pr & (d.q < 0.05)).sum()), pass_bh05_formula=int((pr & (d.q_formula < 0.05)).sum())))
    out = pd.DataFrame(rows); out.to_csv(UT / "08_formel.csv", index=False); print(out.to_string(index=False))
    print("total passing BH 5 %: permutation", out.pass_bh05_permutation.sum(), "formula", out.pass_bh05_formula.sum(), "of", out.hpa_prognostic.sum())

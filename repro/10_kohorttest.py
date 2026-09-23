#!/usr/bin/env python3
"""Whole-cohort permutation test (07) corrected across the 21 cohorts: Bonferroni and Benjamini-Hochberg at 5 %.
p_global_minp from 07 has resolution 1 / (N + 1): 0.005 with 200 permutations, 0.0005 with 2000. Bonferroni needs p < 0.05 / 21 = 0.0024, so a cohort
run with 200 permutations cannot pass it however strong its signal. "unresolved": p is not below the threshold, but the 95 %
Clopper-Pearson interval for the Monte Carlo estimate still reaches below it, so more permutations could change the answer. BH q for a cohort depends only on the
p-values at or above its own, so a coarse p in a stronger cohort does not change q for a weaker one. No randomness.
Output: repro/ut/10_kohorter.csv"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import beta

UT = Path(__file__).resolve().parent / "ut"
s = pd.read_csv(UT / "07_sammanfattning.csv")[["cohort", "events", "permutations", "observed_minp", "null_max", "p_global_minp"]].sort_values("p_global_minp").reset_index(drop=True)
n = len(s); p = s.p_global_minp.to_numpy(); k = np.arange(1, n + 1); bonf = 0.05 / n
s["bh_q"] = np.round(np.minimum.accumulate((p * n / k)[::-1])[::-1], 4); s["bh_pass"] = s.bh_q <= 0.05
hits = np.round(p * (s.permutations + 1)).astype(int) - 1   # permutations at or above the observed count
lower = np.where(hits > 0, beta.ppf(0.025, np.maximum(hits, 1), s.permutations - hits + 1), 0.0)
s["bonferroni"] = np.where(p < bonf, "pass", np.where(lower < bonf, "unresolved", "fail"))
s["observed_above_null_max"] = s.observed_minp > s.null_max
s.to_csv(UT / "10_kohorter.csv", index=False); print(s.to_string(index=False)); print(f"Bonferroni threshold {bonf:.4f}; BH pass: {int(s.bh_pass.sum())} of {n}")

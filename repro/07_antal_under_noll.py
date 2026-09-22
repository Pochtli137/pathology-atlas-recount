#!/usr/bin/env python3
"""How many "prognostic" genes does chance alone produce in ONE dataset? 03 gives the mean (genes x alpha_eff). Genes are correlated, so
the count in a single dataset scatters far more than a binomial would. Here every permutation is run on ALL genes of the cohort, which
keeps the gene-gene correlation intact, and the number of genes below the threshold is recorded per permutation.
  null_count            distribution of #genes with best-cutoff p < 0.001 under permutation (mean, median, 5th, 95th, max)
  p_global              (1 + #permutations with null_count >= observed count) / (N + 1): does the cohort as a whole carry more
                        prognostic genes than chance produces? One test per cohort, valid under any gene-gene dependence.
  the same for univariate Cox at p < 0.001, and the number of genes Cox + Benjamini-Hochberg (q < 0.05) reports under permutation.
Usage: python3 repro/07_antal_under_noll.py [N_PERM]   (default 200 per cohort, run in batches of 10 across processes). Seed fixed.
Output: repro/ut/07_counts_<ABBR>.csv (one row per permutation), repro/ut/07_sammanfattning.csv"""
import sys, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from logrank import prepare, minp_for_genes  # noqa: E402
from cox import cox_for_genes  # noqa: E402
spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
UT = rep.UT; NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 200; BATCH = 10


def bh_count(p, level=0.05):
    n = len(p); ps = np.sort(p); ok = np.flatnonzero(ps <= level * np.arange(1, n + 1) / n); return int(ok[-1] + 1) if ok.size else 0


def batch(arg):
    abbr, b = arg; rng = np.random.default_rng([20260922, sum(map(ord, abbr)), b])
    c, genes, m, tid, handelse, _ = rep.ladda(abbr); expr = m[np.nanmean(m, axis=1) > 1].astype(np.float64); lg = np.log2(np.nan_to_num(expr) + 1); prep = prepare(tid, handelse); rows = []
    for i in range(BATCH):
        perm = rng.permutation(expr.shape[1]); p = minp_for_genes(expr[:, perm], prep)[0]; pc = cox_for_genes(lg[:, perm], prep)["p_lrt"]
        rows.append(dict(cohort=abbr, batch=b, i=i, minp_below_001=int((p < 1e-3).sum()), cox_below_001=int((pc < 1e-3).sum()), cox_bh_q05=bh_count(pc)))
    return rows


if __name__ == "__main__":
    from multiprocessing import Pool
    vikt = {"BRCA": 9, "HNSC": 8, "LUSC": 8, "LUAD": 7, "KIRC": 7, "OV": 6, "STAD": 5, "LIHC": 4}
    tasks = sorted(((a, b) for a in rep.NAMN for b in range(NPERM // BATCH)), key=lambda t: -vikt.get(t[0], 1)); rows = []
    with Pool(11) as pool:
        for k, r in enumerate(pool.imap_unordered(batch, tasks)):
            rows += r
            if k % 20 == 0: print(f"{k + 1}/{len(tasks)} batches", flush=True)
    d = pd.DataFrame(rows).sort_values(["cohort", "batch", "i"]); s5 = pd.read_csv(UT / "05_sammanfattning.csv").set_index("cohort"); out = []
    for a, g in d.groupby("cohort"):
        g.to_csv(UT / f"07_counts_{a}.csv", index=False); x = g.minp_below_001.to_numpy(); xc = g.cox_below_001.to_numpy(); obs, obs_c = int(s5.hpa_prognostic[a]), int(s5.cox_p001[a])
        out.append(dict(cohort=a, events=int(s5.events[a]), permutations=len(g), observed_minp=obs, null_mean=round(x.mean(), 1), null_median=int(np.median(x)), null_p05=int(np.percentile(x, 5)),
                        null_p95=int(np.percentile(x, 95)), null_max=int(x.max()), p_global_minp=round((1 + (x >= obs).sum()) / (len(x) + 1), 4),
                        observed_cox=obs_c, cox_null_mean=round(xc.mean(), 1), cox_null_p95=int(np.percentile(xc, 95)), cox_null_max=int(xc.max()), p_global_cox=round((1 + (xc >= obs_c).sum()) / (len(xc) + 1), 4),
                        cox_bh_any_under_null=round(float((g.cox_bh_q05 > 0).mean()), 3), cox_bh_null_max=int(g.cox_bh_q05.max())))
    out = pd.DataFrame(out); out.to_csv(UT / "07_sammanfattning.csv", index=False); print(out.to_string(index=False))

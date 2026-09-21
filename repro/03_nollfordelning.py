#!/usr/bin/env python3
"""The one extension: what does "p < 0.001" mean under the null when the cutoff is chosen to minimise p?

For each TCGA cohort, patient labels are permuted (breaking any link between expression and survival while keeping both marginals,
including ties in expression and the censoring pattern) and the same best-cutoff procedure is run on a random subset of real genes.
The pooled null min-p values give
  alpha_eff   = P_null(min-p < 0.001)            the real false-positive rate of HPA's threshold in this cohort
  expected_fp = genes tested x alpha_eff         how many "prognostic" genes pure noise would produce
  p_cal       = calibrated p per gene            (1 + #null <= p_obs) / (N + 1)
  q           = Benjamini-Hochberg on p_cal
Usage: python3 repro/03_nollfordelning.py [N_PERM] [GENES_PER_PERM]   (default 40 x 1000 = 40 000 null draws per cohort)
Output: repro/ut/03_null_<ABBR>.npy, repro/ut/03_<ABBR>_gener.csv, repro/ut/03_sammanfattning.csv. Seed fixed."""
import sys, time, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from logrank import prepare, minp_for_genes  # noqa: E402
spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
UT = rep.UT; NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 40; NGEN = int(sys.argv[2]) if len(sys.argv) > 2 else 1000


def bh(p):
    n = len(p); o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1); q = np.minimum.accumulate(q[::-1])[::-1]; ut = np.empty(n); ut[o] = np.minimum(q, 1); return ut


def en_kohort(abbr):
    t0 = time.time(); rng = np.random.default_rng(20260921 + sum(map(ord, abbr)))
    c, genes, m, tid, handelse, _ = rep.ladda(abbr); expr = m[np.nanmean(m, axis=1) > 1].astype(np.float64); prep = prepare(tid, handelse)
    null = []
    for _ in range(NPERM):
        g = rng.choice(expr.shape[0], size=min(NGEN, expr.shape[0]), replace=False); perm = rng.permutation(expr.shape[1])
        null.append(minp_for_genes(expr[g][:, perm], prep)[0])
    null = np.sort(np.concatenate(null)); np.save(UT / f"03_null_{abbr}.npy", null)
    d = pd.read_csv(UT / f"02_{abbr}_gener.csv").dropna(subset=["p_ours"])
    d["p_cal"] = (1 + np.searchsorted(null, d["p_ours"].to_numpy(), side="right")) / (len(null) + 1); d["q"] = bh(d["p_cal"].to_numpy())
    d.to_csv(UT / f"03_{abbr}_gener.csv", index=False)
    a = float((null < 1e-3).mean()); n_prog = int((d["p_ours"] < 1e-3).sum())
    return dict(cohort=abbr, patients=len(tid), events=int(handelse.sum()), genes=len(d), null_draws=len(null), null_median_p=round(float(np.median(null)), 4),
                alpha_eff=round(a, 5), inflation=round(a / 1e-3, 1), prognostic=n_prog, expected_fp=round(len(d) * a, 1), share_expected_fp=round(min(1.0, len(d) * a / max(n_prog, 1)), 3),
                survive_q05=int(((d["p_ours"] < 1e-3) & (d["q"] < 0.05)).sum()), survive_q01=int(((d["p_ours"] < 1e-3) & (d["q"] < 0.01)).sum()), seconds=round(time.time() - t0, 1))


if __name__ == "__main__":
    from multiprocessing import Pool
    with Pool(10) as pool:   # one cohort per process; largest cohorts first
        ordning = sorted(rep.NAMN, key=lambda a: -{"BRCA": 9, "HNSC": 8, "LUSC": 8, "LUAD": 7, "KIRC": 7, "OV": 6, "STAD": 5, "LIHC": 4}.get(a, 1))
        rader = []
        for r in pool.imap_unordered(en_kohort, ordning): rader.append(r); print(r, flush=True)
    pd.DataFrame(rader).to_csv(UT / "03_sammanfattning.csv", index=False)

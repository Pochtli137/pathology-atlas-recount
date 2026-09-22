#!/usr/bin/env python3
"""Is the Cox likelihood-ratio p-value calibrated in the far tail, where Benjamini-Hochberg operates (around 0.05 / 14,000 = 3.6e-6)?
A reviewer suggested the excess BH discoveries under permutation in 07 come from tail miscalibration, not gene-gene correlation.
For each cohort: permute survival, run Cox on all genes, pool the null p-values (about 3 million per cohort), and report the actual
tail rate divided by the nominal threshold. Also records the same for the best-cutoff min-p, corrected with Lausen-Schumacher, since
that is the other p-value fed to BH. Seed fixed. Output: repro/ut/09_cox_svans.csv"""
import sys, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from logrank import prepare
from cox import cox_for_genes
spec = importlib.util.spec_from_file_location("rep", HERE / "02_reproducera.py"); rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 200
T = [1e-3, 1e-4, 1e-5, 3.6e-6, 1e-6]


def en(abbr):
    rng = np.random.default_rng([20260922, 9, sum(map(ord, abbr))]); c, genes, m, tid, ev, _ = rep.ladda(abbr)
    lg = np.log2(np.nan_to_num(m[np.nanmean(m, axis=1) > 1].astype(np.float64)) + 1); prep = prepare(tid, ev); cnt = np.zeros(len(T)); N = 0; minp = 1.0
    for _ in range(NPERM):
        p = cox_for_genes(lg[:, rng.permutation(lg.shape[1])], prep)["p_lrt"]; N += len(p); minp = min(minp, p.min())
        for k, t in enumerate(T): cnt[k] += (p < t).sum()
    return dict(cohort=abbr, events=int(ev.sum()), genes=lg.shape[0], null_draws=N, **{f"ratio_at_{t:g}": round(cnt[k] / N / t, 2) for k, t in enumerate(T)}, min_p=float(minp))


if __name__ == "__main__":
    from multiprocessing import Pool
    with Pool(10) as pool:
        rows = list(pool.imap_unordered(en, sorted(rep.NAMN, key=lambda a: -{"BRCA": 9, "HNSC": 8, "LUSC": 8, "LUAD": 7, "KIRC": 7, "OV": 6}.get(a, 1))))
    d = pd.DataFrame(rows).sort_values("events"); d.to_csv(rep.UT / "09_cox_svans.csv", index=False); print(d.to_string(index=False))

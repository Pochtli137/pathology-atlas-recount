"""Log-rank test for every expression cutoff of many genes at once (the HPA "best cut-off" procedure).

HPA v2 method (Yuan et al., eBioMedicine 2024, Methods): for each gene, every TPM value between the 20th and 80th percentile is tried
as a cutoff that splits patients into a high and a low group, and the cutoff with the lowest log-rank p-value is kept.

The statistic is the standard two-group log-rank chi-square with the hypergeometric variance, as in R's survival::survdiff (rho = 0),
which is what the authors' R code uses. Verified against scipy.stats.logrank in test_logrank.py.

For one cohort with n patients and J distinct event times:
  at_risk[i, j] = 1 if time_i >= t_j            (n x J)
  d_j           = events at t_j, n_j = patients at risk at t_j
Sort patients by expression. The high group for cutoff k is the patients above position k, so its at-risk counts for all cutoffs are a
reverse cumulative sum over the sorted rows. That gives every cutoff's test in O(n * J) per gene instead of O(n * J) per cutoff.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import chi2


def prepare(time: np.ndarray, event: np.ndarray):
    """Cohort-level quantities that do not depend on the gene."""
    time = np.asarray(time, float); event = np.asarray(event, bool)
    tj = np.unique(time[event])                                   # distinct event times
    at_risk = (time[:, None] >= tj[None, :]).astype(np.float64)   # n x J
    died = ((time[:, None] == tj[None, :]) & event[:, None]).astype(np.float64)
    nj = at_risk.sum(0); dj = died.sum(0)
    var_w = np.where(nj > 1, dj * (nj - dj) / np.maximum(nj - 1, 1), 0.0)   # hypergeometric variance weight per event time
    return dict(at_risk=at_risk, died=died, nj=nj, dj=dj, var_w=var_w, n=len(time))


def minp_for_genes(expr: np.ndarray, prep: dict, lo: float = 0.20, hi: float = 0.80, min_group: int = 1):
    """expr: genes x patients. Returns (min_p, best_cutoff_value, n_high_at_best, direction) per gene.
    direction = +1 if the high-expression group has MORE events than expected (unfavourable), -1 if fewer (favourable)."""
    G, n = expr.shape
    A, D, nj, dj, vw = prep["at_risk"], prep["died"], prep["nj"], prep["dj"], prep["var_w"]
    out_p = np.ones(G); out_cut = np.full(G, np.nan); out_nh = np.zeros(G, int); out_dir = np.zeros(G, int)
    for g in range(G):
        x = expr[g]
        order = np.argsort(x, kind="stable"); xs = x[order]
        # high group for cutoff value c = patients with x > c. Cumulate from the top.
        n1 = np.cumsum(A[order][::-1], axis=0)[::-1]      # n1[k, j]: at risk among sorted patients k..n-1
        o1 = np.cumsum(D[order][::-1], axis=0)[::-1]      # observed events among them, per event time
        # candidate cutoffs: distinct values between the 20th and 80th percentile; "high" = strictly greater than the cutoff
        qlo, qhi = np.quantile(x, [lo, hi])
        first_above = np.searchsorted(xs, xs, side="right")            # index of first patient strictly above xs[k]
        cand = np.flatnonzero((xs >= qlo) & (xs <= qhi))
        cand = cand[np.r_[True, xs[cand][1:] != xs[cand][:-1]]] if cand.size else cand   # one per distinct value
        k = first_above[cand]; ok = (k >= min_group) & (n - k >= min_group) & (k < n)
        cand, k = cand[ok], k[ok]
        if k.size == 0: continue
        n1k = n1[k]; O = o1[k].sum(1); E = (n1k * (dj / nj)).sum(1)
        frac = n1k / nj; V = (frac * (1 - frac) * vw).sum(1)
        with np.errstate(divide="ignore", invalid="ignore"): stat = np.where(V > 0, (O - E) ** 2 / V, 0.0)
        b = int(np.argmax(stat))
        out_p[g] = chi2.sf(stat[b], 1); out_cut[g] = xs[cand[b]]; out_nh[g] = n - k[b]; out_dir[g] = 1 if O[b] > E[b] else -1
    return out_p, out_cut, out_nh, out_dir


def logrank_two_groups(time, event, high):
    """Plain two-group log-rank p-value, for testing."""
    prep = prepare(time, event); high = np.asarray(high, bool)
    n1 = prep["at_risk"][high].sum(0); O = prep["died"][high].sum(); E = (n1 * prep["dj"] / prep["nj"]).sum()
    frac = n1 / prep["nj"]; V = (frac * (1 - frac) * prep["var_w"]).sum()
    return float(chi2.sf((O - E) ** 2 / V, 1))

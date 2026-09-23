"""Univariate Cox proportional hazards regression for many genes at once (sensitivity analysis, no cutoff involved).

Model per gene: h(t | x) = h0(t) * exp(beta * x), x = standardised log2(pTPM + 1). Ties are handled with Breslow's approximation.
R's survival::coxph defaults to Efron; with TCGA's day-resolution times the two differ in the third decimal of beta. Breslow is used
because every quantity is then a matrix product with the n x J at-risk matrix from logrank.prepare, so thousands of genes fit at once.

For G genes, W = exp(beta * X) (G x n), A = at-risk (n x J), d_j = events at t_j:
  S0 = W A, S1 = (W*X) A, S2 = (W*X*X) A                       (G x J)
  loglik = beta * sum_events x_i - sum_j d_j log S0_j
  score  = sum_events x_i - sum_j d_j S1_j / S0_j
  info   = sum_j d_j (S2_j / S0_j - (S1_j / S0_j)^2)
Newton-Raphson with step halving. Verified against statsmodels PHReg(ties="breslow") and a brute-force optimiser in test_cox.py.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import chi2


def _parts(beta, X, A, dj, sx):
    W = np.exp(np.clip(beta[:, None] * X, -50, 50)); WX = W * X
    S0 = W @ A; S1 = WX @ A; S2 = (WX * X) @ A
    m = S1 / S0
    ll = beta * sx - (np.log(S0) * dj).sum(1)
    return ll, sx - (m * dj).sum(1), ((S2 / S0 - m * m) * dj).sum(1)


def _parts_strata(beta, X, strata, sx):
    """Stratified partial likelihood: each stratum has its own baseline hazard; loglik, score and information add over strata."""
    ll = beta * sx; u = sx.copy(); info = np.zeros_like(beta)
    for idx, A, dj in strata:
        Xs = X[:, idx]; W = np.exp(np.clip(beta[:, None] * Xs, -50, 50)); WX = W * Xs
        S0 = W @ A; S1 = WX @ A; S2 = (WX * Xs) @ A; m = S1 / S0
        ll = ll - (np.log(S0) * dj).sum(1); u = u - (m * dj).sum(1); info = info + ((S2 / S0 - m * m) * dj).sum(1)
    return ll, u, info


def cox_for_genes(expr: np.ndarray, prep: dict, max_iter: int = 25, tol: float = 1e-9, chunk: int = 2000, strata=None):
    """expr: genes x patients, already on the scale to be modelled (rows with zero variance get p = 1).
    strata: optional list of (patient index array, prep for those patients) for a stratified model (separate baseline hazard per stratum,
    one beta). prep is then only used for the event indicator of all patients.
    Returns dict of arrays per gene: beta (per SD), se, p_lrt, p_wald, p_score, converged."""
    G, n = expr.shape; A, D, dj = prep["at_risk"], prep["died"], prep["dj"]; is_event = D.sum(1) > 0
    if strata is not None:
        st = [(idx, p["at_risk"], p["dj"]) for idx, p in strata]; parts = lambda b, X, sx: _parts_strata(b, X, st, sx)
    else:
        parts = lambda b, X, sx: _parts(b, X, A, dj, sx)
    out = {k: np.full(G, np.nan) for k in ("beta", "se", "p_lrt", "p_wald", "p_score")}; out["converged"] = np.zeros(G, bool)
    for s in range(0, G, chunk):
        X = np.asarray(expr[s:s + chunk], float); sd = X.std(1); ok = sd > 0
        X = (X - X.mean(1, keepdims=True)) / np.where(ok, sd, 1)[:, None]; X[~ok] = 0.0
        sx = X[:, is_event].sum(1); beta = np.zeros(len(X))
        ll0, u0, i0 = parts(beta, X, sx); ll = ll0.copy(); u, info = u0, i0; done = ~ok
        for _ in range(max_iter):
            step = np.where(done | (info <= 0), 0.0, u / np.where(info > 0, info, 1)); new = beta + step
            ll_n, u_n, i_n = parts(new, X, sx)
            for _h in range(10):                                   # step halving where the likelihood went down
                bad = (ll_n < ll - 1e-12) & ~done
                if not bad.any(): break
                step = np.where(bad, step / 2, step); new = beta + step; ll_n, u_n, i_n = parts(new, X, sx)
            done = done | (np.abs(ll_n - ll) < tol * (np.abs(ll) + 1)); beta, ll, u, info = new, ll_n, u_n, i_n
            if done.all(): break
        se = 1 / np.sqrt(np.where(info > 0, info, np.nan)); sl = slice(s, s + len(X))
        out["beta"][sl] = np.where(ok, beta, 0.0); out["se"][sl] = se
        out["p_lrt"][sl] = np.where(ok, chi2.sf(np.maximum(2 * (ll - ll0), 0), 1), 1.0)
        out["p_wald"][sl] = np.where(ok, chi2.sf((beta / se) ** 2, 1), 1.0)
        out["p_score"][sl] = np.where(ok & (i0 > 0), chi2.sf(u0 ** 2 / np.where(i0 > 0, i0, 1), 1), 1.0)
        out["converged"][sl] = done
    return out

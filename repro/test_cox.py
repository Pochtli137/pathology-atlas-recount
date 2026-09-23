"""cox.cox_for_genes against (1) a brute-force optimiser of the Breslow partial likelihood and (2) statsmodels PHReg when installed."""
import numpy as np, pytest
from scipy.optimize import minimize_scalar
from logrank import prepare
from cox import cox_for_genes


def _data(seed=1, n=180, G=12):
    rng = np.random.default_rng(seed); X = rng.normal(size=(G, n)); X[0] = np.round(X[0])          # one gene with heavy ties in x
    haz = np.exp(0.5 * X[1] - 0.3 * X[2]); t = np.ceil(rng.exponential(300 / haz))                  # integer days, so tied event times
    c = np.ceil(rng.exponential(400, n)); return X, np.minimum(t, c), t <= c


def _negll(b, x, time, event):
    ll = 0.0
    for tj in np.unique(time[event]):
        d = ((time == tj) & event); ll += b * x[d].sum() - d.sum() * np.log(np.exp(b * x[time >= tj]).sum())
    return -ll


def test_against_brute_force():
    X, time, event = _data(); r = cox_for_genes(X, prepare(time, event)); assert r["converged"].all()
    for g in range(len(X)):
        x = (X[g] - X[g].mean()) / X[g].std(); b = minimize_scalar(_negll, args=(x, time, event), bounds=(-5, 5), method="bounded", options=dict(xatol=1e-10)).x
        assert abs(b - r["beta"][g]) < 1e-5


def test_against_statsmodels():
    sm = pytest.importorskip("statsmodels.api"); X, time, event = _data(seed=7); r = cox_for_genes(X, prepare(time, event))
    for g in range(len(X)):
        x = (X[g] - X[g].mean()) / X[g].std(); f = sm.PHReg(time, x[:, None], status=event.astype(int), ties="breslow").fit()
        assert abs(f.params[0] - r["beta"][g]) < 1e-5 and abs(f.bse[0] - r["se"][g]) < 1e-5 and abs(f.pvalues[0] - r["p_wald"][g]) < 1e-6


def test_constant_gene_is_p_one():
    X, time, event = _data(); X[3] = 2.0; r = cox_for_genes(X, prepare(time, event)); assert r["p_lrt"][3] == 1.0 and r["beta"][3] == 0.0


def test_stratified_against_statsmodels():
    sm = pytest.importorskip("statsmodels.api"); X, time, event = _data(seed=11, n=240); g_str = np.random.default_rng(3).integers(0, 3, len(time))
    strata = [(np.flatnonzero(g_str == k), prepare(time[g_str == k], event[g_str == k])) for k in range(3)]
    r = cox_for_genes(X, prepare(time, event), strata=strata)
    for g in range(len(X)):
        x = (X[g] - X[g].mean()) / X[g].std(); f = sm.PHReg(time, x[:, None], status=event.astype(int), ties="breslow", strata=g_str).fit()
        assert abs(f.params[0] - r["beta"][g]) < 1e-5 and abs(f.bse[0] - r["se"][g]) < 1e-5


def test_one_stratum_equals_unstratified():
    X, time, event = _data(seed=5); a = cox_for_genes(X, prepare(time, event))
    b = cox_for_genes(X, prepare(time, event), strata=[(np.arange(len(time)), prepare(time, event))])
    assert np.allclose(a["beta"], b["beta"]) and np.allclose(a["p_lrt"], b["p_lrt"])


def test_stratified_against_brute_force():
    X, time, event = _data(seed=11, n=240); g_str = np.random.default_rng(3).integers(0, 3, len(time))
    strata = [(np.flatnonzero(g_str == k), prepare(time[g_str == k], event[g_str == k])) for k in range(3)]
    r = cox_for_genes(X, prepare(time, event), strata=strata); assert r["converged"].all()
    for g in range(len(X)):
        x = (X[g] - X[g].mean()) / X[g].std()
        f = lambda b: sum(_negll(b, x[g_str == k], time[g_str == k], event[g_str == k]) for k in range(3))
        b = minimize_scalar(f, bounds=(-5, 5), method="bounded", options=dict(xatol=1e-10)).x; assert abs(b - r["beta"][g]) < 1e-5

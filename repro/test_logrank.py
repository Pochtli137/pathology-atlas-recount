"""Checks the vectorised log-rank against scipy.stats.logrank and against a brute-force search over cutoffs."""
import numpy as np
from scipy import stats
from logrank import prepare, minp_for_genes, logrank_two_groups


def _data(seed, n=180):
    r = np.random.default_rng(seed)
    t = np.round(r.exponential(1000, n)) + 1; e = r.random(n) < 0.6
    return t, e, r


def test_two_group_matches_scipy():
    for seed in range(5):
        t, e, r = _data(seed); high = r.random(len(t)) < 0.4
        ours = logrank_two_groups(t, e, high)
        a = stats.CensoredData(uncensored=t[high & e], right=t[high & ~e]); b = stats.CensoredData(uncensored=t[~high & e], right=t[~high & ~e])
        ref = stats.logrank(a, b).pvalue
        assert abs(np.log10(ours) - np.log10(ref)) < 1e-6, (seed, ours, ref)


def test_minp_matches_brute_force():
    t, e, r = _data(11); x = r.lognormal(2, 1, (4, len(t))); prep = prepare(t, e)
    p, cut, nh, _ = minp_for_genes(x, prep)
    for g in range(4):
        lo, hi = np.quantile(x[g], [0.2, 0.8]); best = 1.0
        for c in np.unique(x[g][(x[g] >= lo) & (x[g] <= hi)]):
            high = x[g] > c
            if high.sum() and (~high).sum(): best = min(best, logrank_two_groups(t, e, high))
        assert abs(np.log10(p[g]) - np.log10(best)) < 1e-6, (g, p[g], best)


def test_minp_is_anticonservative_under_null():
    """The point of the whole project in one test: with no signal at all, min-p over cutoffs rejects far more often than its nominal level."""
    t, e, r = _data(3, n=300); x = r.normal(size=(400, 300)); p, *_ = minp_for_genes(x, prepare(t, e))
    assert (p < 0.05).mean() > 0.25 and np.median(p) < 0.2

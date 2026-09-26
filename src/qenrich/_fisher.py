"""Vectorized hypergeometric p-values for ORA, matching clusterProfiler's test.

``clusterProfiler::enrichGO`` runs ``phyper(k - 1, M, N - M, n, lower.tail =
FALSE)`` for every gene set that contains at least one query gene: a one-sided
over-representation test, P(X >= k) with X ~ Hypergeom(N, M, n). This module
computes the same numbers for every term of a gene set in one pass.

Why not call scipy per table: in the builds measured here ``hypergeom.pmf``/
``cdf`` cost 40-150 us per call (``betaln`` costs 0.06 us per element), and the
per-table route costs ~0.5 ms per (term, set), i.e. ~30s for a three-set run on a
30k-gene net. That was the whole runtime.

The pmf goes through ``betaln`` and the CDF sums a recurrence anchored at the
hypergeometric mode. Agreement with R's ``phyper`` and ``enrichit::ora_gson`` is
~1e-10 relative on a 30k-gene net, ~3e-12 on small ones; the residual is
special-function rounding, not an algorithmic difference.
"""

import numpy as np
from scipy.special import betaln


def hg_logpmf(k, K, n, N):
    """log P(X = k) for X ~ Hypergeom(N, K, n); the formula scipy's _logpmf uses.

    ``betaln`` evaluates in log space, so nothing underflows until the value
    itself is denormal.
    """
    bad = N - K
    return (betaln(K + 1, 1) + betaln(bad + 1, 1) + betaln(N - n + 1, n + 1)
            - betaln(k + 1, K - k + 1) - betaln(n - k + 1, bad - n + k + 1)
            - betaln(N + 1, 1))


def hg_cdf(k, K, n, N):
    """P(X <= k), walking outward from the mode so no partial sum underflows.

    Starting at the lower support end instead computes ``exp(logpmf(lo))``, which
    is exactly 0.0 once lo sits far in the tail (logpmf can reach -600), and every
    later term is then 0 * ratio: the CDF comes out 0.0 where the true value is 1.0.
    The mode pmf is bounded below by roughly 1/(sigma*sqrt(2*pi)) with
    sigma <= sqrt(N/4), so for any N that fits in int64 the anchor stays positive.

    k below the support gives 0, at or above the top gives 1, matching scipy.
    """
    k = np.asarray(k, np.int64)
    lo = np.maximum(0, n - (N - K))
    hi = np.minimum(K, n)
    mode = np.clip((n + 1) * (K + 1) // (N + 2), lo, hi)
    target = np.clip(k, lo, hi)
    anchor = np.exp(hg_logpmf(mode, K, n, N))
    cum = np.where(target >= mode, anchor, 0.0)

    # down from the mode: p(j) = p(j+1) * (j+1)(N-K-n+j+1) / ((K-j)(n-j))
    j, p = mode - 1, anchor.copy()
    live = j >= lo
    while live.any():
        i = np.where(live)[0]
        jj, Kk, nn, NN = j[i], K[i], n[i], N[i]
        pj = p[i] * ((jj + 1) * (NN - Kk - nn + jj + 1)) / ((Kk - jj) * (nn - jj))
        cum[i] += np.where(jj <= target[i], pj, 0.0)
        p[i], j[i] = pj, jj - 1
        live[i] = j[i] >= lo[i]

    # up from the mode: p(j) = p(j-1) * (K-j+1)(n-j+1) / (j(N-K-n+j)); terms past
    # the target are never needed, so those tests drop out early
    j, p = mode + 1, anchor.copy()
    live = (j <= hi) & (j <= target)
    while live.any():
        i = np.where(live)[0]
        jj, Kk, nn, NN = j[i], K[i], n[i], N[i]
        pj = p[i] * ((Kk - jj + 1) * (nn - jj + 1)) / (jj * (NN - Kk - nn + jj))
        cum[i] += pj
        p[i], j[i] = pj, jj + 1
        live[i] = (j[i] <= hi[i]) & (j[i] <= target[i])

    cum = np.clip(cum, 0.0, 1.0)
    cum[k < lo] = 0.0
    cum[k >= hi] = 1.0
    return cum


def fisher_pvalues(K, n, N, obs, alternative="greater"):
    """Hypergeometric p-values for many tables sharing the same margins.

    Parameters
    ----------
    K, n, N : array_like of int
        Term size, set size and universe size per test (broadcast against ``obs``).
    obs : array_like of int
        Observed overlap per test, i.e. the ``a`` cell of ``[[a, b], [c, d]]``.
    alternative : {'greater', 'less'}
        ``greater`` (default) is P(X >= k), the over-representation test
        clusterProfiler runs. ``less`` is P(X <= k), which tests depletion.

    Returns
    -------
    np.ndarray of float
        One p-value per test.
    """
    K = np.asarray(K, dtype=np.int64)
    n = np.asarray(n, dtype=np.int64)
    N = np.asarray(N, dtype=np.int64)
    obs = np.asarray(obs, dtype=np.int64)

    if alternative == "greater":
        # R: phyper(k - 1, M, N - M, n, lower.tail = FALSE). Written as the CDF of
        # the complementary hypergeometric, which is how scipy computes the same
        # upper tail and keeps the floats close to both references.
        return hg_cdf(K - obs, K, N - n, N)
    if alternative == "less":
        return hg_cdf(obs, K, n, N)
    raise ValueError(f"alternative must be greater or less: {alternative!r}")

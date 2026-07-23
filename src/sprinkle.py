"""
Shared utilities: 2D (1+1-dimensional) Minkowski Poisson sprinkling, and the
causal-order-nearest-predecessor statistic R_i = u_i = dx/dt used throughout
Section 6 of the paper as the concrete, checkable equivariant reconstruction
ingredient.

Coordinates: (t, x), metric ds^2 = -dt^2 + dx^2. Causal order: j precedes i
(j -> i) iff t_i - t_j > 0 and |x_i - x_j| < t_i - t_j. The "nearest" causal
predecessor of i is the j preceding i with minimal proper-time interval
d_tau^2 = dt^2 - dx^2 (both dt and d_tau^2 are Lorentz-invariant quantities,
which is exactly why this statistic is causal-order-equivariant in the sense
of Section 6.1 of the paper).
"""

import numpy as np


def sprinkle(T: float, X: float, density: float, rng: np.random.Generator) -> np.ndarray:
    """
    Poisson-sprinkle points uniformly in the box t in [-T,T], x in [-X,X],
    with the given number-density (points per unit spacetime area). Since
    dt*dx is the Lorentz-invariant Minkowski volume form in 1+1D, a uniform
    Poisson process in these coordinates is exactly the standard
    Poincare-invariant-in-law causal-set sprinkling (Bombelli, Lee, Meyer &
    Sorkin 1987).
    """
    n = rng.poisson(density * (2 * T) * (2 * X))
    t = rng.uniform(-T, T, n)
    x = rng.uniform(-X, X, n)
    return np.stack([t, x], axis=1)


def lorentz_boost(pts: np.ndarray, beta: float) -> np.ndarray:
    """Apply an exact Lorentz boost with velocity parameter beta (|beta|<1)."""
    gamma = 1.0 / np.sqrt(1.0 - beta ** 2)
    t, x = pts[:, 0], pts[:, 1]
    t2 = gamma * (t + beta * x)
    x2 = gamma * (x + beta * t)
    return np.stack([t2, x2], axis=1)


def nearest_causal_predecessor_bruteforce(pts: np.ndarray, query_idx=None):
    """
    O(n * |query_idx|) brute-force version: for each queried point i, search
    ALL other points for the causal predecessor of minimal proper-time
    interval. Returns (u, pred_index) arrays aligned with query_idx.
    Suitable for exact small-scale checks (few thousand points); use
    nearest_causal_predecessor_fast for larger, windowed searches.
    """
    t, x = pts[:, 0], pts[:, 1]
    idxs = query_idx if query_idx is not None else np.arange(len(pts))
    u = np.full(len(idxs), np.nan)
    pred = np.full(len(idxs), -1, dtype=int)
    for k, i in enumerate(idxs):
        dt = t[i] - t
        dx = x[i] - x
        causal = (dt > 1e-9) & (np.abs(dx) < dt)
        if not np.any(causal):
            continue
        dtau2 = dt[causal] ** 2 - dx[causal] ** 2
        j_local = np.argmin(dtau2)
        candidate_idx = np.where(causal)[0][j_local]
        pred[k] = candidate_idx
        u[k] = dx[causal][j_local] / dt[causal][j_local]
    return u, pred


def nearest_causal_predecessor_fast(pts: np.ndarray, query_idx: np.ndarray, window: float):
    """
    Efficient version for large point sets: restricts the predecessor search
    to points within `window` of the query point's t-coordinate (via a
    sorted-array slice), rather than scanning all n points. `window` must be
    chosen much larger than the typical predecessor scale (~1/sqrt(density))
    and much larger than any boost-induced shear of the sampling region, or
    results will be biased (see boundary_artifact_demo.py / README).
    """
    t, x = pts[:, 0], pts[:, 1]
    order = np.argsort(t)
    t_sorted, x_sorted = t[order], x[order]

    u = np.full(len(query_idx), np.nan)
    for k, i in enumerate(query_idx):
        ti, xi = t[i], x[i]
        lo = np.searchsorted(t_sorted, ti - window, side="left")
        hi = np.searchsorted(t_sorted, ti, side="left")
        if hi <= lo:
            continue
        tc, xc = t_sorted[lo:hi], x_sorted[lo:hi]
        dt = ti - tc
        dx = xi - xc
        causal = (dt > 1e-9) & (np.abs(dx) < dt)
        if not np.any(causal):
            continue
        dtau2 = dt[causal] ** 2 - dx[causal] ** 2
        j = np.argmin(dtau2)
        u[k] = dx[causal][j] / dt[causal][j]
    return u


def relativistic_velocity_addition(u: np.ndarray, beta: float) -> np.ndarray:
    """u' = (u + beta) / (1 + u*beta), the exact special-relativistic
    aberration/velocity-addition formula that an equivariant local velocity
    statistic must satisfy under a boost."""
    return (u + beta) / (1.0 + u * beta)

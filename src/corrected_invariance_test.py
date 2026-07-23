"""
Corrected version of the distributional boost-invariance test (Section 6.3,
Table 2 rows 4-5): restricts query points to an interior region of the
sampling box, with margin much larger than both the boost-induced shear and
the typical causal-predecessor length scale (~1/sqrt(density)), so that
censoring artifacts (see boundary_artifact_demo.py) cannot occur.

Compares:
  (a) native-vs-boosted: KS test between the unboosted statistic distribution
      and the boosted-frame statistic distribution from an independent draw.
  (b) native-vs-native baseline: KS test between two independent unboosted
      draws, as the reference for "statistically indistinguishable."

If Lemma D (Section 6.1) is correct, (a) and (b) should give comparable KS
statistics and p-values -- i.e. the boosted-frame distribution should be no
more distinguishable from native than two independent native draws are from
each other.

Uses the windowed (efficient) predecessor search since this test requires a
much larger point set than the O(n^2) exact check.

Usage:
    python3 corrected_invariance_test.py --seed 7
"""

import argparse
import numpy as np
from scipy import stats

from sprinkle import sprinkle, lorentz_boost, nearest_causal_predecessor_fast


def draw_and_measure(T, X, density, margin, window, n_query, rng, beta=None):
    """Draw a fresh sprinkling, optionally boost it, and compute u on a set
    of query points restricted to the interior (margin away from all edges
    of the UNBOOSTED box, which is what matters: the search window is small
    compared to margin, so the predecessor search never needs points outside
    the sampled box, whether or not a boost is subsequently applied)."""
    pts = sprinkle(T, X, density, rng)
    core = np.where((np.abs(pts[:, 0]) < T - margin) & (np.abs(pts[:, 1]) < X - margin))[0]
    q = rng.choice(core, size=min(n_query, len(core)), replace=False)

    if beta is not None:
        pts_used = lorentz_boost(pts, beta)
    else:
        pts_used = pts

    u = nearest_causal_predecessor_fast(pts_used, q, window=window)
    return u[~np.isnan(u)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--X", type=float, default=300.0)
    ap.add_argument("--density", type=float, default=3.0)
    ap.add_argument("--margin", type=float, default=40.0,
                     help="min distance of query points from box edges (>> window)")
    ap.add_argument("--window", type=float, default=20.0,
                     help="predecessor search window (>> 1/sqrt(density))")
    ap.add_argument("--n_query", type=int, default=4000)
    ap.add_argument("--beta", type=float, default=0.6)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    uA = draw_and_measure(args.T, args.X, args.density, args.margin, args.window,
                           args.n_query, rng, beta=None)
    uB_boosted = draw_and_measure(args.T, args.X, args.density, args.margin, args.window,
                                   args.n_query, rng, beta=args.beta)
    uC = draw_and_measure(args.T, args.X, args.density, args.margin, args.window,
                           args.n_query, rng, beta=None)

    ks_main = stats.ks_2samp(uA, uB_boosted)
    ks_base = stats.ks_2samp(uA, uC)

    print("CORRECTED TEST (interior margin control)")
    print("-" * 70)
    print(f"n(A, native) = {len(uA)}   n(B, boosted) = {len(uB_boosted)}   n(C, native) = {len(uC)}")
    print(f"mean(uA)          = {uA.mean():.4f}   std(uA)          = {uA.std():.4f}")
    print(f"mean(uB, boosted) = {uB_boosted.mean():.4f}   std(uB, boosted) = {uB_boosted.std():.4f}")
    print()
    print(f"KS (native vs. boosted) : D={ks_main.statistic:.4f}, p={ks_main.pvalue:.4f}")
    print(f"KS (native vs. native)  : D={ks_base.statistic:.4f}, p={ks_base.pvalue:.4f}   <- baseline")
    print()
    if ks_main.pvalue > 0.05 and abs(ks_main.statistic - ks_base.statistic) < 0.03:
        print("Result: boosted-frame distribution is statistically indistinguishable")
        print("from the native distribution, consistent with Lemma D.")
    else:
        print("Result: distributions differ beyond the native-vs-native baseline --")
        print("check margin/window parameters before concluding a real effect.")

    return dict(ks_main=ks_main._asdict() if hasattr(ks_main, "_asdict") else
                dict(statistic=float(ks_main.statistic), pvalue=float(ks_main.pvalue)),
                ks_base=dict(statistic=float(ks_base.statistic), pvalue=float(ks_base.pvalue)))


if __name__ == "__main__":
    main()

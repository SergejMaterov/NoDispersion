"""
Reproduces the SPURIOUS violation of boost-invariance found in a naive
distributional test (Section 6.3, Table 2 row 3): sampling query points from
a fixed rectangular box without any margin, then boosting and comparing
distributions, appears to show a violation (KS test rejects, p ~ 0) because a
Lorentz boost maps the rectangular sampling box into a sheared parallelogram,
censoring some true causal predecessors near the (now slanted) edges.

The box here is deliberately sized so that essentially every point is
"close to an edge" relative to the box size -- this is precisely the
regime that triggers the artifact; see corrected_invariance_test.py for
the same check with a properly controlled margin.

This script is deliberately kept as a companion to
corrected_invariance_test.py: running both side by side is the point --
it demonstrates the artifact and its resolution, rather than only the
resolution.

Usage:
    python3 boundary_artifact_demo.py --seed 7
"""

import argparse
import numpy as np
from scipy import stats

from sprinkle import sprinkle, lorentz_boost, nearest_causal_predecessor_bruteforce


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=20.0, help="small box -> no margin from edges (deliberate)")
    ap.add_argument("--X", type=float, default=20.0)
    ap.add_argument("--density", type=float, default=2.0)
    ap.add_argument("--beta", type=float, default=0.6)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    # Sample A: native, unboosted -- ALL points in the (small) box act as both
    # the population and the query set (no margin control at all).
    ptsA = sprinkle(args.T, args.X, args.density, rng)
    uA, _ = nearest_causal_predecessor_bruteforce(ptsA)
    uA = uA[~np.isnan(uA)]

    # Sample B: fresh draw, boosted -- same small-box, no-margin setup.
    ptsB = sprinkle(args.T, args.X, args.density, rng)
    ptsB_boosted = lorentz_boost(ptsB, args.beta)
    uB, _ = nearest_causal_predecessor_bruteforce(ptsB_boosted)
    uB = uB[~np.isnan(uB)]

    ks_stat, ks_p = stats.ks_2samp(uA, uB)

    print("NAIVE TEST (no margin control) -- expect a SPURIOUS apparent violation")
    print("-" * 70)
    print(f"n(A, native) = {len(uA)}   n(B, boosted) = {len(uB)}")
    print(f"mean(u, native)         = {uA.mean():.4f}")
    print(f"mean(u, boosted, naive) = {uB.mean():.4f}   <-- shifted away from 0: artifact, not physics")
    print(f"KS two-sample statistic : D={ks_stat:.4f}, p={ks_p:.4g}")
    print()
    print("See corrected_invariance_test.py for the resolution: restricting")
    print("query points to a region well inside the sampling box removes this")
    print("effect entirely.")

    return dict(mean_native=float(uA.mean()), mean_boosted_naive=float(uB.mean()),
                ks_stat=float(ks_stat), ks_p=float(ks_p))


if __name__ == "__main__":
    main()

"""
Exact (non-statistical) verification that the causal-order-nearest-predecessor
statistic is causal-order-equivariant (Section 6.1-6.2 of the paper, Table 2
rows 1-2):

  1. The identity of the nearest causal predecessor is preserved EXACTLY under
     a Lorentz boost (because causal order and proper-time interval are
     themselves Lorentz-invariant quantities).
  2. The transformed local velocity u' matches the exact relativistic
     velocity-addition formula u' = (u+beta)/(1+u*beta) to machine precision.

Usage:
    python3 exact_equivariance_check.py --n 3000 --beta 0.6 --seed 42
"""

import argparse
import numpy as np

from sprinkle import (sprinkle, lorentz_boost, nearest_causal_predecessor_bruteforce,
                       relativistic_velocity_addition)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T", type=float, default=400.0)
    ap.add_argument("--X", type=float, default=400.0)
    ap.add_argument("--density", type=float, default=2.0)
    ap.add_argument("--n", type=int, default=3000, help="number of points to subsample for the O(n^2) brute-force check")
    ap.add_argument("--beta", type=float, default=0.6)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    pts = sprinkle(args.T, args.X, args.density, rng)
    sub = rng.choice(len(pts), size=min(args.n, len(pts)), replace=False)
    pts_sub = pts[sub]

    u_orig, pred_orig = nearest_causal_predecessor_bruteforce(pts_sub)

    pts_boosted = lorentz_boost(pts_sub, args.beta)
    u_boost, pred_boost = nearest_causal_predecessor_bruteforce(pts_boosted)

    frac_same_pred = np.nanmean((pred_orig == pred_boost).astype(float))

    valid = (~np.isnan(u_orig)) & (~np.isnan(u_boost))
    u_predicted = relativistic_velocity_addition(u_orig[valid], args.beta)
    max_err = float(np.nanmax(np.abs(u_predicted - u_boost[valid])))

    print(f"n points checked        : {valid.sum()}")
    print(f"beta                     : {args.beta}")
    print(f"predecessor identity preserved (fraction): {frac_same_pred:.6f}  (expect 1.0 exactly)")
    print(f"max |u_boost - relativistic_velocity_addition(u_orig)| : {max_err:.3e}  (expect ~machine precision)")

    return dict(n=int(valid.sum()), beta=args.beta, frac_same_pred=frac_same_pred, max_err=max_err)


if __name__ == "__main__":
    main()

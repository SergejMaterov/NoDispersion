"""
Sanity-check (not a proof -- Chebyshev's inequality is proved analytically in
the paper, Section 6.5) that the finite-N bound

    P( |Phi_N| > k*sqrt(N*v) ) <= 1/k^2

holds empirically for a simple uncorrelated-phase-contribution model, and
that it is loose (as Chebyshev bounds generically are) compared to the true
tail probability once a specific distribution is assumed -- illustrating why
the bound is stated as a distribution-free guarantee rather than a tight
estimate.

Usage:
    python3 chebyshev_sanity_check.py --N 10000 --trials 200000
"""

import argparse
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000, help="number of traversed configurations")
    ap.add_argument("--trials", type=int, default=200_000, help="number of independent Phi_N realizations")
    ap.add_argument("--v", type=float, default=1.0, help="Var(delta phi)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    # delta_phi_i drawn i.i.d. mean-zero, variance v (uniform, so NOT Gaussian --
    # illustrates the bound holds beyond the Gaussian case, as Chebyshev requires
    # only finite variance, not any specific distribution or full independence).
    half_width = np.sqrt(3 * args.v)  # uniform[-hw,hw] has variance v
    deltas = rng.uniform(-half_width, half_width, size=(args.trials, args.N))
    Phi_N = deltas.sum(axis=1)

    sigma_N = np.sqrt(args.N * args.v)

    print(f"N={args.N}, trials={args.trials}, v={args.v}")
    print(f"empirical std(Phi_N) = {Phi_N.std():.4f}   theoretical sqrt(N*v) = {sigma_N:.4f}")
    print()
    print(f"{'k':>4}  {'Chebyshev bound 1/k^2':>22}  {'empirical P(|Phi_N|>k*sigma)':>28}")
    for k in (1, 2, 3, 5, 10):
        bound = 1.0 / k ** 2
        empirical = float(np.mean(np.abs(Phi_N) > k * sigma_N))
        print(f"{k:>4}  {bound:>22.4%}  {empirical:>28.4%}")

    print()
    print("Expectation: empirical probability <= Chebyshev bound at every k")
    print("(the bound is distribution-free and therefore conservative/loose;")
    print("this is expected and does not indicate an error).")


if __name__ == "__main__":
    main()

# Computational Appendix — Boost Invariance of Decorrelation (Section 6)

This repository reproduces every numerical claim in **Section 6 ("Boost Invariance of Decorrelation,
Strengthened: A Rigorous Sub-Case and an Explicit Conjecture")** of:

> S. Materov, *On the Absence of Automatic Linear Lorentz-Violating Dispersion from Discreteness
> Alone*, 2026.

It supersedes the original Theorem 2 ("proof sketch", not independently checked) with:

1. A **rigorous lemma** (Lemma D) for the sub-class of models built on a Poincaré-invariant-in-law
   point process (in particular, causal sets with Lorentz-invariant Poisson sprinkling), proved from
   two independently-established facts rather than a circular appeal to "no preferred frame."
2. **Exact algebraic verification** that a concrete, causal-order-equivariant local statistic
   transforms exactly as required under a boost (machine-precision agreement).
3. A **documented methodological pitfall**: a naive numerical test of boost-invariance produces a
   spurious violation from finite-box boundary censoring, and the corrected test that resolves it.
4. A **quantitative Chebyshev bound** added to Theorem 1, with an empirical sanity check.

## Repository structure

```
.
├── README.md
├── requirements.txt
├── run_all.sh                          <- reproduces every row of Table 2 + the Theorem 1 addendum
├── LICENSE.md
├── src/
│   ├── sprinkle.py                     <- shared point-process / causal-predecessor utilities
│   ├── exact_equivariance_check.py     <- Table 2, rows 1-2 (exact algebraic check)
│   ├── boundary_artifact_demo.py       <- Table 2, row 3 (the spurious violation, deliberately reproduced)
│   ├── corrected_invariance_test.py    <- Table 2, rows 4-5 (the fix + baseline)
│   └── chebyshev_sanity_check.py       <- Section 6.5 addendum to Theorem 1
└── results/
    └── table2_log.txt                  <- full output log from run_all.sh
```

## Why the "naive" script is deliberately kept in the repository

`boundary_artifact_demo.py` reproduces a **wrong-looking result on purpose**. Early testing found an
apparent boost-invariance violation; the cause (not a real effect) was a rectangular sampling box
becoming a sheared parallelogram under a boost, which censors some true causal predecessors near the
now-slanted edges. This is kept as a runnable script, not just a paragraph in the paper, because it is
a **generic trap for anyone else attempting to numerically test (EM4) or its boost-invariance** in a
specific discrete-gravity model (CDT, spinfoam, etc.) — running it and seeing the spurious violation
first-hand is more convincing (and more useful as a warning) than being told about it.

`corrected_invariance_test.py` is the same test with the fix: query points are restricted to a region
well inside the sampling box, with margin much larger than both the boost-induced shear and the typical
causal-predecessor length scale. Run both, in that order, to see the artifact appear and then disappear.

## Quickstart

```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
bash run_all.sh
```

Runtime: well under a minute in total. Individual scripts can also be run directly with `python3
src/<script>.py --help` to see tunable parameters (box size, density, boost velocity, margin, seed).

## Key results (Table 2 in the paper)

| Check | Method | Result |
|---|---|---|
| Predecessor identity preserved under boost | Exact algebraic (2957 points, β=0.6) | 100.0% identical |
| u_boost = relativistic velocity addition of u_orig | Exact algebraic | max error ~4×10⁻¹³ (machine precision) |
| Naive test, small box, no margin control | KS two-sample, native vs. boosted | D≈0.10, p≈0 (spurious — artifact) |
| Corrected test, large box + interior margin | KS two-sample, native vs. boosted | D≈0.02, p≈0.30 |
| Baseline: native vs. independent native draw | KS two-sample | D≈0.02, p≈0.47 |

Minor run-to-run variation (a few percent in D, more in p) is expected from the random seed and
sample size and does not affect the qualitative conclusion: the corrected KS statistic is
indistinguishable from the native-vs-native baseline, while the naive one is not.

## What each script proves, and what it doesn't

- `exact_equivariance_check.py` is a **deterministic** check (given the seed): it proves the chosen
  statistic transforms exactly as the relativistic velocity-addition formula predicts, for the actual
  points drawn. It says nothing about the *distribution* of the statistic — that is the job of the
  other two scripts.
- `boundary_artifact_demo.py` and `corrected_invariance_test.py` are statistical (Monte Carlo + KS
  test) and therefore have sampling noise; the qualitative conclusion (naive test fails, corrected test
  matches baseline) is robust across seeds, but exact D/p values will vary slightly run to run.
- `chebyshev_sanity_check.py` illustrates, but does not by itself prove, the analytic Chebyshev bound
  in Section 6.5; the proof is a two-line application of Chebyshev's inequality to Var(Φ_N)=Nv, given
  in the paper. The script exists to show the bound is respected (and generically loose, as
  distribution-free bounds are) for a concrete, non-Gaussian example.

## Reproducibility notes

- All scripts take an explicit `--seed`; the values used in the paper are `seed=42`
  (`exact_equivariance_check.py`) and `seed=7` (`boundary_artifact_demo.py`,
  `corrected_invariance_test.py`).
- `corrected_invariance_test.py` uses a windowed (sorted-array) predecessor search
  (`nearest_causal_predecessor_fast` in `sprinkle.py`) so it can run efficiently on point sets of
  ~10⁶ points; `exact_equivariance_check.py` and `boundary_artifact_demo.py` use the brute-force
  O(n²) search (`nearest_causal_predecessor_bruteforce`) since they intentionally work with smaller
  point sets (a few thousand points) where an independent, unoptimized implementation is preferable
  for a "no shortcuts" sanity check.

## Citing

If this computational appendix is cited independently of the paper, please cite the paper itself
(DOI: [10.5281/zenodo.21915193](https://doi.org/10.5281/zenodo.21915193) and reference this repository as its computational supplement.

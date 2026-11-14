# Computational Supplement — Conditional GOE Benchmark

This repository contains the code and numerical data supporting Sections 6.2–6.3 of:

> S. Materov, **Symmetry Constraints on Quantitative Intra-Multiplet Splittings** (2026).

The calculations are a **conditional statistical benchmark**. They assume that a real-symmetric perturbation matrix on a real `d`-dimensional multiplet is assigned an `O(d)`-invariant maximum-entropy prior at fixed Frobenius second moment. This additional orthogonal isotropy is not implied by invariance under an arbitrary physical symmetry group.

## Scope and normalization

The matrix distribution is the Gaussian Orthogonal Ensemble (GOE) with independent upper-triangular entries

- diagonal entries: `N(0, σ²)`;
- off-diagonal entries: `N(0, σ²/2)`.

The statistic studied is the global absolute-eigenvalue ratio

`R = max_i |λ_i| / min_i |λ_i|`.

It is a statistic of the perturbation matrix. It is **not automatically a physical particle-mass ratio**; that interpretation requires an additional model connecting the perturbation matrix to the mass or mass-squared operator.

## Results

### Exact result for `d = 2`

For `GOE(2)`, the exact tail distribution is

`P(R > R₀) = sqrt(2 / (1 + R₀²))`, for `R₀ ≥ 1`.

Consequently:

| Quantity | Exact value |
|---|---:|
| Median | `sqrt(7) ≈ 2.6458` |
| 90th percentile | `sqrt(199) ≈ 14.1067` |
| 99th percentile | `sqrt(19999) ≈ 141.4178` |
| `P(R > 20)` | `sqrt(2/401) ≈ 7.0622%` |
| `P(R > 200)` | `sqrt(2/40001) ≈ 0.7071%` |

The exact result is checked in two independent ways: direct GOE(2) Monte Carlo sampling and deterministic numerical quadrature over the joint `(m, r)` density. The quadrature verification concerns the **two-dimensional exact result**; it is not a `d = 3` Gauss–Hermite calculation.

### Monte Carlo benchmark for `d = 3, 4, 5, 6`

For dimensions above two, the repository reports Monte Carlo estimates, not exact quantiles. The scan uses 3,000,000 matrices per dimension and seed `12345`. Bootstrap standard errors are reported for the sample medians (200 subsampled bootstrap resamples of size at most 200,000).

| `d` | Median (bootstrap SE) | 90th pct. | 99th pct. | `P(R > 20)` | `P(R > 200)` |
|---:|---:|---:|---:|---:|---:|
| 2 (exact) | 2.646 | 14.107 | 141.42 | 7.06% | 0.707% |
| 3 | 5.006 (0.013) | 26.87 | 271.47 | 13.42% | 1.35% |
| 4 | 7.226 (0.020) | 39.10 | 392.30 | 19.36% | 1.96% |
| 5 | 9.561 (0.025) | 51.73 | 522.65 | 25.44% | 2.60% |
| 6 | 11.814 (0.029) | 63.77 | 636.69 | 31.11% | 3.19% |

The increasing median is an extreme-value effect: the statistic takes the maximum and minimum absolute eigenvalues over a growing number of eigenvalues. These values characterize this GOE prior only; they do not define a model-independent typical mass hierarchy and do not follow from symmetry invariance alone.

## Repository structure

```text
.
├── README.md
├── requirements.txt
├── run_all.sh
├── LICENSE.md
├── src/
│   ├── exact_formulas.py
│   ├── quadrature_verification.py
│   ├── monte_carlo_d2.py
│   └── general_d_scan.py
└── results/
    └── general_d_table.csv
```

## Reproduction

Create an environment and install the dependencies:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
bash run_all.sh
```

The main individual scripts are:

```bash
python3 src/exact_formulas.py
python3 src/quadrature_verification.py
python3 src/monte_carlo_d2.py --n 20000000 --seed 0
python3 src/general_d_scan.py --n 3000000 --seed 12345 --out results/general_d_table.csv
```

`general_d_scan.py` reports the sample median, its subsampled-bootstrap standard error, sample 90th and 99th percentiles, and empirical tail probabilities. The table in `results/general_d_table.csv` is the data source for the manuscript's `d = 3` through `d = 6` rows.

## Citing

Please cite the manuscript and identify this repository as its computational supplement:

S. Materov, *Symmetry Constraints on Quantitative Intra-Multiplet Splittings* (2026), computational supplement, <https://github.com/SergejMaterov/SymmetryNoGo>.

#!/usr/bin/env bash
# Reproduces every numerical claim in Section 6 ("Boost Invariance of
# Decorrelation, Strengthened") of the paper, in order. Run from the
# repository root. Total runtime: well under a minute.
set -euo pipefail

OUT=results/table2_log.txt
mkdir -p results
: > "$OUT"

echo "############################################################" | tee -a "$OUT"
echo "# 1-2. Exact equivariance check (Table 2, rows 1-2)"          | tee -a "$OUT"
echo "############################################################" | tee -a "$OUT"
python3 src/exact_equivariance_check.py --n 3000 --beta 0.6 --seed 42 | tee -a "$OUT"

echo | tee -a "$OUT"
echo "############################################################" | tee -a "$OUT"
echo "# 3. Naive test: spurious violation from box-boundary artifact (row 3)" | tee -a "$OUT"
echo "############################################################" | tee -a "$OUT"
python3 src/boundary_artifact_demo.py --seed 7 | tee -a "$OUT"

echo | tee -a "$OUT"
echo "############################################################" | tee -a "$OUT"
echo "# 4-5. Corrected test with margin control + baseline (rows 4-5)" | tee -a "$OUT"
echo "############################################################" | tee -a "$OUT"
python3 src/corrected_invariance_test.py --seed 7 | tee -a "$OUT"

echo | tee -a "$OUT"
echo "############################################################" | tee -a "$OUT"
echo "# 6. Chebyshev bound sanity check (Section 6.5 addendum to Theorem 1)" | tee -a "$OUT"
echo "############################################################" | tee -a "$OUT"
python3 src/chebyshev_sanity_check.py --N 5000 --trials 50000 | tee -a "$OUT"

echo | tee -a "$OUT"
echo "Done. Full log written to $OUT"

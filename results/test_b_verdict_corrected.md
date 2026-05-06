# Test B verdict (corrected criterion 3)

## Why this re-evaluation exists

The original criterion 3 required `autocorrelation_peak_lag == N` for every cell. This is mathematically wrong: a sequence Y with period N is also periodic with period 2N, 3N, ..., so the autocorrelation function has equal-magnitude peaks at every integer multiple of N. Sampling noise determines which multiple wins the argmax. The framework's actual prediction is *Y is periodic with period N* (i.e., y_{i+N} = y_i in the noiseless case).

Corrected operationalization:

- **Test 3a:** `peak_lag % N == 0` (peak at some integer multiple of N)
- **Test 3b:** `|peak_val - predicted_match_prob| < 0.05`, with `predicted_match_prob = (1-eps)^2 + eps^2/(N-1)`. When peak_lag is a multiple of N, peak_val equals the lag-N autocorrelation up to sampling noise, so the peak_val column is a valid stand-in.
- **Skip:** cells with `L == N` (only one polymer cycle, lag-N autocorr has zero pairs).

## Summary of all three criteria

- **Criterion 1** (low-noise saturation, I_empirical -> log_2(N) at eps=0.001): PASS (established by prior analysis; every cell within 0.00-0.01% relative error)
- **Criterion 2** (no L-growth of saturated information): PASS (established by prior analysis; every (N, eps) group has max/min ratio < 1.0001)
- **Criterion 3 (corrected)** (peak at multiple of N AND peak_val matches predicted match prob within 0.05): PASS

## Criterion 3 breakdown

- Total cells: 105
- Skipped (L == N, no second cycle): 21
- Eligible cells (L > N): 84
- Test 3a passed (peak_lag is multiple of N): 84 / 84
- Test 3b passed (peak_val within 0.05 of predicted): 84 / 84
- Both tests passed: 84 / 84

## Remaining criterion-3 failures

None. All eligible cells pass both 3a and 3b.

## Skipped cells (L == N)

These have only one polymer cycle, so lag-N autocorrelation has zero
pairs to average. The autocorrelation engine returned a small-lag peak
with low value (consistent with no detectable periodicity from a single cycle).

| N | eps | L | peak_lag | peak_val |
|---|-----|---|----------|----------|
| 2 | 0.001 | 2 | 1 | 0.00195 |
| 2 | 0.01 | 2 | 1 | 0.02012 |
| 2 | 0.05 | 2 | 1 | 0.09577 |
| 3 | 0.001 | 3 | 2 | 0.00107 |
| 3 | 0.01 | 3 | 1 | 0.00983 |
| 3 | 0.05 | 3 | 1 | 0.04886 |
| 4 | 0.001 | 4 | 3 | 0.00078 |
| 4 | 0.01 | 4 | 2 | 0.00664 |
| 4 | 0.05 | 4 | 2 | 0.03262 |
| 5 | 0.001 | 5 | 4 | 0.00052 |
| 5 | 0.01 | 5 | 2 | 0.00491 |
| 5 | 0.05 | 5 | 4 | 0.02466 |
| 6 | 0.001 | 6 | 4 | 0.00040 |
| 6 | 0.01 | 6 | 2 | 0.00407 |
| 6 | 0.05 | 6 | 5 | 0.01963 |
| 8 | 0.001 | 8 | 1 | 0.00029 |
| 8 | 0.01 | 8 | 6 | 0.00293 |
| 8 | 0.05 | 8 | 7 | 0.01402 |
| 10 | 0.001 | 10 | 2 | 0.00023 |
| 10 | 0.01 | 10 | 2 | 0.00223 |
| 10 | 0.05 | 10 | 7 | 0.01094 |

## Final verdict

# **Test B: OVERALL PASS**

All three criteria pass under the corrected criterion-3 logic. The framework's prediction that Y is periodic with period N is borne out: peaks land at integer multiples of N (3a), and the lag-N autocorrelation matches the noise-model prediction `(1-eps)^2 + eps^2/(N-1)` to within sampling tolerance (3b).

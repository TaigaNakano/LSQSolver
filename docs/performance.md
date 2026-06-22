# Performance

This page summarizes preliminary benchmark results for LSQSolver and GNU Octave.

The benchmark uses dense square matrices with sizes

`n = 10, 50, 100, 500, 1000, 2000`

and two rank conditions:

* `rank_half`: square matrices with numerical rank `n / 2`
* `full_rank`: square matrices with numerical rank `n`

Each value below is the median of 10 runs.
All timings are reported in milliseconds.

## Benchmark Environment

The benchmark was run on the following machine:

| Item             | Value                                               |
| ---------------- | --------------------------------------------------- |
| Model            | MacBook Pro                                         |
| Chip             | Apple M1 Pro                                        |
| CPU Cores        | 8 total: 6 performance cores and 2 efficiency cores |
| Memory           | 16 GB                                               |

## LSQSolver

| case      |      10 |      50 |     100 |      500 |     1000 |      2000 |
| --------- | ------: | ------: | ------: | -------: | -------: | --------: |
| rank_half | 10.4537 | 3.15045 | 6.41995 |   97.064 | 380.2049 | 1701.5698 |
| full_rank |  0.1438 |  0.3684 | 1.05525 | 23.15535 | 153.6233 | 1134.0057 |

## GNU Octave

The following table shows Octave solve timings for `backslash`, `pinv`, and QR factorization.

| case / method         |  10 |  50 |  100 |    500 |    1000 |    2000 |
| --------------------- | --: | --: | ---: | -----: | ------: | ------: |
| rank_half / backslash | 0.5 | 0.1 |  0.2 |    6.2 |   27.15 |  132.05 |
| rank_half / pinv      | 0.1 | 0.6 | 3.95 | 204.55 |  1240.9 | 9750.85 |
| rank_half / qr_factor | 0.1 | 0.6 |  2.7 |   80.9 |   483.2 |  2102.2 |
| full_rank / backslash | 0.0 | 0.1 |  0.2 |   5.15 |   23.65 |   106.2 |
| full_rank / pinv      | 0.1 | 0.8 | 4.75 |  322.6 | 1788.95 | 15599.7 |
| full_rank / qr_factor | 0.0 | 0.5 |  2.2 |   67.8 |   258.1 | 1482.95 |

## Comparison with Octave QR Factorization

Since LSQSolver is based on column-pivoted QR decomposition, the closest Octave comparison is `qr_factor`.

| case / solver                |      10 |      50 |     100 |      500 |     1000 |      2000 |
| ---------------------------- | ------: | ------: | ------: | -------: | -------: | --------: |
| rank_half / LSQSolver        | 10.4537 | 3.15045 | 6.41995 |   97.064 | 380.2049 | 1701.5698 |
| rank_half / Octave qr_factor |     0.1 |     0.6 |     2.7 |     80.9 |    483.2 |    2102.2 |
| full_rank / LSQSolver        |  0.1438 |  0.3684 | 1.05525 | 23.15535 | 153.6233 | 1134.0057 |
| full_rank / Octave qr_factor |     0.0 |     0.5 |     2.2 |     67.8 |    258.1 |   1482.95 |

## Observations

For very small matrices, overhead dominates the timing results. 

For full-rank matrices, LSQSolver is faster than Octave QR factorization for `n >= 50` in this benchmark.

For rank-deficient matrices, LSQSolver is slower than Octave QR factorization at smaller sizes, but becomes competitive for larger matrices. In this benchmark, LSQSolver is faster than Octave QR factorization for `n = 1000` and `n = 2000`.

Compared with Octave `pinv`, LSQSolver is significantly faster for large matrices. This is expected because `pinv` is typically based on singular value decomposition, while LSQSolver avoids SVD and uses QR-based rank detection with Cholesky-based minimum-norm reconstruction.

## Notes

These results are preliminary and depend on hardware, runtime, compiler settings, BLAS/LAPACK configuration, and benchmark implementation details.

The Octave `qr_factor` timing measures QR factorization behavior and does not necessarily include the full minimum-norm solution reconstruction performed by LSQSolver for rank-deficient cases. Therefore, the comparison should be interpreted as a practical reference rather than a strict one-to-one algorithmic comparison.

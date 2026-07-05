# Performance

This page summarizes preliminary benchmark results for LSQSolver and GNU Octave.

The benchmark uses dense square matrices with sizes

`n = 10, 50, 100, 500, 1000, 2000`

and two rank conditions:

- `rank_half`: square matrices with numerical rank `n / 2`
- `full_rank`: square matrices with numerical rank `n`

Each value below is the median of 10 runs.
All timings are reported in milliseconds.

For the `rank_half` 10 x 10 case, the shorter post-warm-up timing was used because the first execution was dominated by benchmark start-up overhead.

## Benchmark Environment

The benchmark was run on the following machine:

| Item      | Value                                               |
| --------- | --------------------------------------------------- |
| Model     | MacBook Pro                                         |
| Chip      | Apple M1 Pro                                        |
| CPU Cores | 8 total: 6 performance cores and 2 efficiency cores |
| Memory    | 16 GB                                               |

## LSQSolver

| case      |      10 |    50 |   100 |   500 |  1000 |   2000 |
| --------- | ------: | ----: | ----: | ----: | ----: | -----: |
| rank_half |   0.130 |  1.81 |  6.25 |  57.9 | 213.7 |  764.7 |
| full_rank |  0.0859 | 0.312 | 0.853 |  17.1 | 103.1 |  564.4 |

## GNU Octave

The following table shows Octave timings for `backslash`, `pinv`, and QR factorization.

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

| case / solver                |     10 |    50 |   100 |   500 |  1000 |   2000 |
| ---------------------------- | -----: | ----: | ----: | ----: | ----: | -----: |
| rank_half / LSQSolver        |  0.130 |  1.81 |  6.25 |  57.9 | 213.7 |  764.7 |
| rank_half / Octave qr_factor |    0.1 |   0.6 |   2.7 |  80.9 | 483.2 | 2102.2 |
| full_rank / LSQSolver        | 0.0859 | 0.312 | 0.853 |  17.1 | 103.1 |  564.4 |
| full_rank / Octave qr_factor |    0.0 |   0.5 |   2.2 |  67.8 | 258.1 | 1483.0 |

The following speedups are computed as

`Octave qr_factor time / LSQSolver time`.

Values larger than 1 mean that LSQSolver was faster.

| case      |   10 |   50 |  100 |  500 | 1000 | 2000 |
| --------- | ---: | ---: | ---: | ---: | ---: | ---: |
| rank_half | 0.77 | 0.33 | 0.43 | 1.40 | 2.26 | 2.75 |
| full_rank |    - | 1.60 | 2.58 | 3.97 | 2.50 | 2.63 |

The `full_rank` 10 x 10 speedup is omitted because Octave reported `0.0 ms` for that case.

## Comparison with Octave pinv

| case / solver          |     10 |    50 |   100 |    500 |   1000 |    2000 |
| ---------------------- | -----: | ----: | ----: | -----: | -----: | ------: |
| rank_half / LSQSolver  |  0.130 |  1.81 |  6.25 |   57.9 |  213.7 |   764.7 |
| rank_half / Octave pinv|    0.1 |   0.6 |  3.95 | 204.55 | 1240.9 | 9750.85 |
| full_rank / LSQSolver  | 0.0859 | 0.312 | 0.853 |   17.1 |  103.1 |   564.4 |
| full_rank / Octave pinv|    0.1 |   0.8 |  4.75 |  322.6 | 1789.0 | 15599.7 |

The following speedups are computed as

`Octave pinv time / LSQSolver time`.

| case      |   10 |   50 |  100 |  500 |  1000 |  2000 |
| --------- | ---: | ---: | ---: | ---: | ----: | ----: |
| rank_half | 0.77 | 0.33 | 0.63 | 3.54 |  5.81 | 12.75 |
| full_rank | 1.16 | 2.56 | 5.57 | 18.9 | 17.35 | 27.64 |

## Comparison Across LSQSolver Versions

The following tables compare the elapsed times across several LSQSolver implementations.

### Rank-deficient cases

| version  |      10 |      50 |    100 |     500 |    1000 |     2000 |
| -------- | :------ | :-----: | :----: | :-----: | :-----: | :------: |
| 1.0.0, 1.0.1 | - | 3.15045 | 6.41995 | 97.0640 | 380.205 | 1701.570 |
| 1.0.3 (1.0.2) | - | 3.54720 | 6.98700 | 96.9781 | 338.588 | 1353.362 |
| latest   |  0.130 | 1.80985 | 6.24700 | 57.8635 | 213.684 |  764.698 |

### Full-rank cases

| version  |      10 |      50 |     100 |     500 |    1000 |     2000 |
| :------: | :-----: | :-----: | :-----: | :-----: | :-----: | :------: |
| 1.0.0, 1.0.1 | 0.14380 | 0.36840 | 1.05525 | 23.1554 | 153.623 | 1134.006 |
|1.0.3(1.0.2) | 0.08455 | 0.41210 | 1.09120 | 23.0219 | 148.378 | 1171.723 |
| latest   | 0.08590 | 0.31245 | 0.85320 | 17.0662 | 103.104 |  564.353 |

## Speedup from the Initial Version

The following speedups are computed as

`initial time / latest time`.

Values larger than 1 mean that the latest version was faster.

| case      |    10 |   50 |  100 |  500 | 1000 | 2000 |
| --------- | ----: | ---: | ---: | ---: | ---: | ---: |
| rank_half | 80.38 | 1.74 | 1.03 | 1.68 | 1.78 | 2.23 |
| full_rank |  1.67 | 1.18 | 1.24 | 1.36 | 1.49 | 2.01 |

The `rank_half` 10 x 10 speedup is not representative because the initial measurement included benchmark start-up overhead.

## Large-Case Summary

| case                |  1.0.0, 1.0.1 | 1.0.3 (1.0.2) | latest | initial / latest |
| ------------------- | -----: | ------------: | -----: | ---------------: |
| rank_half, n = 500  |   97.1 |          97.0 |   57.9 |            1.68x |
| rank_half, n = 1000 |  380.2 |         338.6 |  213.7 |            1.78x |
| rank_half, n = 2000 | 1701.6 |        1353.4 |  764.7 |            2.23x |
| full_rank, n = 500  |   23.2 |          23.0 |   17.1 |            1.36x |
| full_rank, n = 1000 |  153.6 |         148.4 |  103.1 |            1.49x |
| full_rank, n = 2000 | 1134.0 |        1171.7 |  564.4 |            2.01x |

## Observations

For very small matrices, benchmark start-up overhead and runtime warm-up can dominate the measured time. In particular, the first `rank_half` 10 x 10 execution showed a large overhead, while the post-warm-up timing was around 0.13 ms.

For large dense matrices, LSQSolver was substantially faster than Octave QR factorization in this benchmark. At `n = 2000`, LSQSolver was about 2.75x faster for the rank-deficient case and about 2.63x faster for the full-rank case.

Compared with Octave `pinv`, LSQSolver was significantly faster for large matrices. At `n = 2000`, LSQSolver was about 12.75x faster for the rank-deficient case and about 27.64x faster for the full-rank case.

Across LSQSolver versions, the latest implementation improved the `n = 2000` elapsed time by about 2.23x for the rank-deficient case and about 2.01x for the full-rank case compared with the initial implementation.

## Notes

These results are preliminary and depend on hardware, runtime, compiler settings, BLAS/LAPACK configuration, and benchmark implementation details.

The Octave `qr_factor` timing measures QR factorization behavior and does not necessarily include the full minimum-norm solution reconstruction performed by LSQSolver for rank-deficient cases. Therefore, the comparison should be interpreted as a practical reference rather than a strict one-to-one algorithmic comparison.

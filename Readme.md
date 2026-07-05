# LSQSolver

[<u>日本語版</u>](https://github.com/TaigaNakano/LSQSolver/blob/main/Readme_jp.md)

A lightweight least-squares solver for dense matrices in .NET.

Supports:

- Overdetermined systems
- Underdetermined systems
- Rank-deficient systems
- Minimum-norm solutions

without SVD and without external dependencies.

GitHub repository: <https://github.com/TaigaNakano/LSQSolver>

---

## Features

- Column-Pivoted QR Decomposition (CPQR)
- Automatic numerical rank detection
- Minimum-norm solutions
- Rank-deficient problem support
- Parallelized factorization using `Parallel.For`
- Optional low-allocation overwrite mode
- Optional storage of QR intermediates (`R`, `Qᵀb`, and pivot information)
- No external dependencies

---
## Why LSQSolver?

Many numerical libraries solve least-squares problems through SVD-based pseudo-inverses.

While SVD is extremely robust, it can be more expensive than necessary when the goal is simply to obtain a least-squares solution.

LSQSolver is based on:

1. Column-Pivoted QR Decomposition (CPQR)
2. Numerical rank detection
3. Cholesky-based minimum-norm reconstruction

This approach supports rank-deficient and underdetermined problems without requiring a full singular value decomposition.

---

## Installation

```bash
dotnet add package LSQSolver
```

## Usage

### Import module

```csharp
using LSQSolver;
using static LSQSolver.LSQSolver;
```

### Construct a matrix

```csharp
double[][] data =
{
    new double[] { 1, 2 },
    new double[] { 3, 4 },
    new double[] { 5, 6 }
};

var A = new MatrixObject(data);
```

### Solve a least-squares problem

```csharp
double[] b = { 7, 8, 9 };
var result = Solve(A, b); // or LSQSolver.LSQSolver.Solve(A, b);
```

### Solve Method

```csharp
var result = Solve(
    A,
    b,
    overwrite: true,
    store_intermediates: false,
    rank_tolerance: 2.22044604925032e-16,
    check_finite: true);
```

### Parameters

- `A`: Coefficient matrix.
- `b`: Right-hand side vector.
- `overwrite`: If `true`, `A` and `b` are overwritten to reduce memory allocation.
- `store_intermediates`: If `true`, stores `R`, `Qᵀb`, and pivot information.
- `rank_tolerance`: Relative tolerance used for numerical rank detection.
- `check_finite`: If `true`, checks whether `A` and `b` contain `NaN` or `Infinity`.

### Return value: LSQSolverResult

`Solve()` returns an `LSQSolverResult` object.

The result object contains the computed solution, solver status, basic diagnostics, and optional QR intermediate data.

Before using `Solution`, check `Status`.

| Property | Type | Description |
| --- | --- | --- |
| `Status` | `LSQSolverStatus` | Status code of the solve operation. Check this value before using `Solution`. |
| `Solution` | `double[]` | Computed least-squares solution vector. For rank-deficient or underdetermined problems, this is the minimum-norm solution when the solve succeeds. |
| `Rows` | `int` | Number of rows of the input matrix `A`. |
| `Cols` | `int` | Number of columns of the input matrix `A`. |
| `Rank` | `int` | Estimated numerical rank detected during the column-pivoted QR factorization. |
| `ResidualNorm` | `double` | Euclidean norm of the residual, `\|\|Ax - b\|\|`. |
| `R` | `MatrixObject?` | Upper-triangular factor obtained by column-pivoted QR factorization. Columns are stored in pivoted order. This is stored only when `store_intermediates` is `true`; otherwise it is `null`. |
| `Qtb` | `double[]?` | Transformed right-hand side vector `Qᵀb`. This is stored only when `store_intermediates` is `true`; otherwise it is `null`. |
| `Pivot` | `int[]?` | Column pivot information. `Pivot[j]` gives the original column index of the `j`-th pivoted column. This is stored only when `store_intermediates` is `true`; otherwise it is `null`. |
| `Tag` | `object?` | Optional tag field reserved for additional metadata. |

#### Public methods

| Method | Description |
| --- | --- |
| `ToString(bool omit = false)` | Converts the result object to a readable string. If `omit` is `true`, long arrays and matrices are shortened for display. |

#### Example

```csharp
var result = Solve(A, b, store_intermediates: true);

if (result.Status != LSQSolverStatus.Success)
{
    Console.WriteLine(result.Status);
    return;
}

double[] x = result.Solution;

Console.WriteLine($"Rank: {result.Rank}");
Console.WriteLine($"Residual norm: {result.ResidualNorm}");
Console.WriteLine(result.ToString(omit: true));
```

## Documentation

The following documents are planned and will be added in future releases.

### Theory

[`theory.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/theory.md)

Topics:

- Least-squares formulation
- Rank-revealing QR decomposition
- Numerical rank detection
- Minimum-norm solutions
- Cholesky-based reconstruction

### Examples

[`polynomial-fit.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/polynomial-fit.md)

- Polynomial curve fitting
- Practical fitting examples

[`gravity-inversion.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/gravity-inversion.md)

- Gravity anomaly inversion
- Underdetermined least-squares problems
- Minimum-norm reconstruction

### Performance

Preliminary benchmark results are available in:

[`performance.md`](https://github.com/TaigaNakano/LSQSolver/blob/main/docs/performance.md)

The benchmark compares LSQSolver with GNU Octave for dense square least-squares problems under both full-rank and rank-deficient conditions.

The reported timings are medians of 10 runs and were measured on the following machine:

| Item | Value |
|---|---|
| Model | MacBook Pro |
| Chip | Apple M1 Pro |
| CPU Cores | 8 total: 6 performance cores and 2 efficiency cores |
| Memory | 16 GB |

In this benchmark, LSQSolver showed the following results:

- For full-rank dense matrices, LSQSolver was faster than Octave QR factorization for `n >= 50` in this benchmark. 
- At `n = 2000`, LSQSolver took `564.4 ms` for the full-rank case, while Octave QR factorization took `1483.0 ms`. 
- For rank-deficient matrices, LSQSolver became faster than Octave QR factorization for larger sizes. At `n = 2000`, LSQSolver took `764.7 ms`, while Octave QR factorization took `2102.2 ms`. 
- Compared with Octave `pinv`, LSQSolver was significantly faster for large matrices. At `n = 2000`, LSQSolver was about `12.8x` faster for the rank-deficient case and about `27.6x` faster for the full-rank case.

These results are preliminary and depend on hardware, runtime, compiler settings, BLAS/LAPACK configuration, and benchmark implementation details.

## Future Plans

Future development will depend on user needs.

Possible extensions include support for regularized least-squares problems, multiple right-hand sides, and weighted least-squares problems. These features may be considered when they can reduce user-side implementation errors, avoid unnecessary memory allocation, or provide algorithmic advantages inside the solver.

## License

MIT License

# LSQSolver

A lightweight least-squares solver for dense matrices in .NET.

Supports:

- Overdetermined systems
- Underdetermined systems
- Rank-deficient systems
- Minimum-norm solutions

without SVD and without external dependencies.

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
var result = LSQSolver.Solve(A, b);
```

### Solve Method

```csharp
var result = LSQSolver.Solve(
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

### Return Value

`Solve()` returns an `LSQSolverResult`.

Check `result.Status` before using `result.Solution`.

```csharp
if (result.Status == LSQSolverStatus.Success)
{
    double[] x = result.Solution;
}
```

### Check solver status

```csharp
if (result.Status != LSQSolverStatus.Success)
{
    Console.WriteLine(result.Status);
    return;
}
```

#### Status Codes

Possible status values include:

- `Success`
- `NullMatrix`
- `EmptyMatrix`
- `NullVector`
- `DimensionMismatch`
- `InvalidMatrix`
- `InvalidVector`
- `CholeskyFailed`


### Access the solution and diagnostics

```csharp
double[] x = result.Solution;

Console.WriteLine(result.Rank);
Console.WriteLine(result.ResidualNorm);
```

## Documentation

The following documents are planned and will be added in future releases.

### Theory

`docs/theory.md`

Topics:

- Least-squares formulation
- Rank-revealing QR decomposition
- Numerical rank detection
- Minimum-norm solutions
- Cholesky-based reconstruction

### Examples

`docs/polynomial-fit.md`

- Polynomial curve fitting
- Vandermonde matrices
- Practical fitting examples

`docs/gravity-inversion.md`

- Gravity anomaly inversion
- Underdetermined least-squares problems
- Minimum-norm reconstruction

### Performance

Benchmark comparisons are planned and will be published in a future release.

Planned comparisons include:

- GNU Octave
- Various matrix sizes and rank conditions
- Execution time

Status: Work in progress.

## License

MIT License

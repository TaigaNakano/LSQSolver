# Changelog

All notable changes to LSQSolver are documented in this file.

The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version numbers follow a simple semantic versioning style.

## [1.0.4] - 2026-07-05

### Added

- Added `CHANGELOG.md` to document version history in the repository.
- Added a gravity inversion example document.
  - The example explains gravity inversion as an underdetermined least-squares problem.
  - It demonstrates how LSQSolver returns a minimum-norm solution and why this is a modeling choice rather than a universal geological prior.
  - Added English and Japanese documentation and generated example figures for shallow and deep anomaly cases.
- Added and expanded theory documentation.
  - The theory page now summarizes the least-squares problem classes handled by LSQSolver.
  - It also documents the overall algorithm flow, including CPQR, rank detection, triangular solve, and minimum-norm reconstruction.
- Updated performance documentation.
  - Added latest benchmark results for rank-deficient and full-rank dense square problems.
  - Added comparisons against GNU Octave QR factorization and `pinv`.
  - Added cross-version performance comparisons and speedup summaries.
- Added generated plot source and figure assets for documentation examples.

### Changed

- Optimized several core numerical kernels in `LSQSolver.cs`.
  - Introduced `System.Runtime.InteropServices` and low-level `Unsafe` / `MemoryMarshal` access patterns in hot loops.
  - Reworked Householder trailing-column updates to process columns in blocks of four inside `Parallel.For`, with scalar handling for the remaining tail columns.
  - Applied similar low-level access patterns to column-norm updates, Householder application, back substitution, and related internal loops.
- Updated README and documentation links to reflect the expanded documentation set.
- Updated polynomial fitting documentation and associated plot assets.
- Updated the solution file structure for the current public repository layout.

### Performance

- In the latest benchmark, LSQSolver took approximately `764.7 ms` for the `rank_half`, `n = 2000` case and `564.4 ms` for the `full_rank`, `n = 2000` case.
- Compared with the initial implementation, the latest version improved the `n = 2000` elapsed time by about `2.23x` for the rank-deficient case and about `2.01x` for the full-rank case.
- Compared with Octave `pinv`, LSQSolver was about `12.75x` faster for the rank-deficient `n = 2000` case and about `27.64x` faster for the full-rank `n = 2000` case in the reported benchmark.

## [1.0.3] - 2026-06-28

### Changed

- Updated the package version from `1.0.2` to `1.0.3`.
- Kept the public solver API unchanged from the preceding package update.

## [1.0.2] - 2026-06-28

### Changed

- Changed the default value of `overwrite` in `Solve()` from `true` to `false`.
  - This makes the default behavior safer because the input matrix `A` and right-hand side vector `b` are not overwritten unless explicitly requested.
- Improved the minimum-norm reconstruction path for rank-deficient and underdetermined problems.
  - The construction of the small symmetric positive definite systems used in the Cholesky-based reconstruction was parallelized.
- Cleaned up the solution structure for the public package.
  - Test and benchmark projects were removed from the public solution file, leaving the library project as the main package target.
- Removed unnecessary diagnostic output from `LSQSolverResult.ToString()`.

## [1.0.1] - 2026-06-21

### Changed

- Updated the README based on `1.0.0`.
- Improved the project documentation while keeping the initial solver implementation unchanged.

## [1.0.0] - 2026-06-21

### Added

- Initial public version of LSQSolver.
- Added a lightweight dense least-squares solver for .NET.
- Added support for:
  - Overdetermined systems
  - Underdetermined systems
  - Rank-deficient systems
  - Minimum-norm solutions
- Added a solver based on column-pivoted QR decomposition.
- Added automatic numerical rank detection.
- Added Cholesky-based minimum-norm reconstruction for rank-deficient and underdetermined cases.
- Added optional overwrite mode to reduce memory allocation.
- Added optional storage of QR intermediate data:
  - `R`
  - `Qᵀb`
  - Pivot information
- Added `MatrixObject` as the matrix container used by the solver.
- Added `LSQSolverResult` for returning the solution, status, diagnostics, and optional intermediate data.
- Added `LSQSolverStatus` for non-throwing status reporting.

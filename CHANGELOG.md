# Changelog

All notable changes to LSQSolver are documented in this file.

The format is inspired by [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version numbers follow a simple semantic versioning style.

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


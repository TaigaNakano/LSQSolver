//using System;
//using System.Threading.Tasks;
using System.Runtime.CompilerServices;

namespace LSQSolver
{
    public static class LSQSolver
    {
        /// <summary>
        /// Unit relative rounding error 
        /// </summary>
        const double EPS = 2.22044604925032e-16;

        /// <summary>
        /// Solves the least squares problem Ax ≈ b using column-pivoted QR decomposition.
        /// Supports full-rank and rank-deficient cases with minimal norm solution.
        /// </summary>
        /// <param name="A">The input matrix A.</param>
        /// <param name="b">The right-hand side vector b.</param>
        /// <param name="overwrite">
        /// If true, A and b are overwritten in place to reduce allocation and memory traffic.
        /// In this mode, the original A and b must not be reused after calling Solve.
        /// If false, A and b are copied before factorization.
        /// </param>
        /// <param name="store_intermediates">
        /// If true, stores QR intermediates useful for follow-up computations: R, Q^T b, and pivot.
        /// </param>
        /// <returns>The computed least squares result. Check Status before using Solution.</returns>
        public static LSQSolverResult Solve(MatrixObject A, double[] b, bool overwrite = false, bool store_intermediates = false, double rank_tolerance =  EPS, bool check_finite = true)
        {
            // Design note:
            // This method intentionally keeps the main numerical workflow in one place.
            // The goal is to reduce allocation, abstraction overhead, and function-call overhead
            // in the hot path. Low-level private kernels may trade readability for speed.
            //
            // Algorithm outline:
            //   1. Column-pivoted Householder QR with GEQP3-ish norm update.
            //   2. Rank detection by pivot column norm.
            //   3. Solve the leading triangular system R11 x1 = Q^T b.
            //   4. If rank-deficient or under-determined, compute the minimum 2-norm completion.
            //   5. Unpivot and return result, optional intermediates, and diagnostics.

            // ------------------------------------------------------------
            // 1. Validate input and initialize workspace
            // ------------------------------------------------------------
            LSQSolverStatus status = ValidateSolveInput(A, b, check_finite);
            if (status != LSQSolverStatus.Success)
            {
                return new LSQSolverResult
                {
                    Status = status,
                    Rows = A?.Rows ?? 0,
                    Cols = A?.Cols ?? 0
                };
            }

            int rows = A.Rows;
            int cols = A.Cols;

            // overwrite == true is the fast path: A and b are intentionally destroyed.
            MatrixObject R = overwrite ? A : new MatrixObject(A);
            double[] Qtb = overwrite ? b : (double[])b.Clone();

            int[] ipiv = new int[cols];
            for (int j = 0; j < cols; j++) ipiv[j] = j;

            double[] x = new double[cols];
            double[] flatten_matrix = R.array; // internal, same assembly

            LSQSolverResult result = new()
            {
                Status = LSQSolverStatus.Success,
                Rows = rows,
                Cols = cols
            };

            // ------------------------------------------------------------
            // 2. Rank-revealing column-pivoted QR
            // ------------------------------------------------------------
            double[] vn1 = new double[cols];
            double[] vn2 = new double[cols];
            InitializeColumnNorms(flatten_matrix, rows, cols, ipiv, vn1, vn2);

            double init_max_norm_of_cols = 0.0;
            for (int j = 0; j < cols; j++)
                if (vn1[j] > init_max_norm_of_cols) init_max_norm_of_cols = vn1[j];

            // Rank criterion based on a NumPy/SciPy-style tolerance scale:
            //
            //   tol = eps * max(m, n) * scale
            //
            // SVD-based routines usually use the largest singular value as scale.
            // This QR-based solver uses the initial maximum column norm instead.
            double rank_tol = rank_tolerance * Math.Max(rows, cols) * init_max_norm_of_cols;
            rank_tol = Math.Max(rank_tol, EPS); // enforce a sanity lower bound to avoid pathological cases
            CPQR(flatten_matrix, rows, cols, ipiv, Qtb, vn1, vn2, rank_tol, out int rankR);
            result.Rank = rankR;    

            // Store QR intermediates before minimum-norm completion.
            // CompleteMinimumNormSolution overwrites the R12 block with Y = R11^{-1}R12.
            if (store_intermediates)
            {
                result.R = new MatrixObject(R);
                result.Qtb = (double[])Qtb.Clone();
                result.Pivot = (int[])ipiv.Clone();
            }

            // ------------------------------------------------------------
            // 3. Backward substitution: x1 = R11^{-1} f1
            // ------------------------------------------------------------
            
            int[] base_rank = new int[rankR];
            for (int j = 0; j < rankR; j++)
                base_rank[j] = ipiv[j] * rows;
            
            BackwardSubstitutionForRHS(flatten_matrix, base_rank, Qtb, x);

            // Compute residual before minimum-norm completion.
            // The completion step may overwrite R12 with Y = R11^{-1}R12.
            double residualNorm = GetResidualNorm(flatten_matrix, rows, cols, Qtb, ipiv, x); // x is still in pivot order
            result.ResidualNorm = residualNorm;

            // ------------------------------------------------------------
            // 4. Rank-deficient / under-determined: minimal 2-norm completion
            // ------------------------------------------------------------
            int available_rank = rows < cols ? rows : cols;
            if (rankR < available_rank || rows < cols)
            {
                if (!CompleteMinimumNormSolution(flatten_matrix, rows, cols, ipiv, base_rank, x))
                    result.Status = LSQSolverStatus.CholeskyFailed;  
            }

            // ------------------------------------------------------------
            // 5. Unpivot and finalize result
            // ------------------------------------------------------------
            result.Solution = UnpivotSolution(cols, ipiv, x);
           
            return result;
        }


        /// <summary>
        /// Validates public Solve inputs without throwing exceptions.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static LSQSolverStatus ValidateSolveInput(MatrixObject? A, double[]? b, bool check_finite = false)
        {
            if (A is null) return LSQSolverStatus.NullMatrix;
            if (A.Rows <= 0 || A.Cols <= 0) return LSQSolverStatus.EmptyMatrix;
            if (b is null) return LSQSolverStatus.NullVector;
            if (b.Length != A.Rows) return LSQSolverStatus.DimensionMismatch;
            if (check_finite)
            {
                for (int i = 0; i < b.Length; i++)
                {
                    if (!double.IsFinite(b[i]))
                        return LSQSolverStatus.InvalidVector;
                }

                double[] arr = A.array;
                for (int i = 0; i < arr.Length; i++)
                {
                    if (!double.IsFinite(arr[i]))
                        return LSQSolverStatus.InvalidMatrix;
                }
            }
            return LSQSolverStatus.Success;
        }

#region  QR private kernels
        /// <summary>
        /// GEQP3-ish norms
        /// </summary>
        /// <param name="A"></param>
        /// <param name="ipiv"></param>
        /// <param name="vn1"></param>
        /// <param name="vn2"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void InitializeColumnNorms(double[] valarr, int rows, int cols,  int[] ipiv, double[] vn1, double[] vn2)
        {
       

            Parallel.For(0, cols, pos =>
            {
                int col = ipiv[pos];
                int basec = col * rows;

                double sum = 0.0;
                for (int i = 0; i < rows; i++)
                {
                    double e = valarr[basec + i];
                    sum += e * e;
                }

                double norm = Math.Sqrt(sum);
                vn1[pos] = norm;
                vn2[pos] = norm;
            });
        }

        private static void CPQR(double[] flatten_matrix, int rows, int cols, int[] ipiv, double[] Qtb, double[] vn1, double[] vn2, double rank_tol, out int rankR)
        {
            int available_rank = rows < cols ? rows : cols;
            rankR = 0;
            for (int j = 0; j < available_rank; j++)
            {
                if (PivotAndCheckRankDeficient(ipiv, vn1, vn2, j, rank_tol)) break;
                ApplyHouseholderToColumn(flatten_matrix, rows, cols, ipiv, j, Qtb);
                UpdateTrailingColumnNorms(flatten_matrix, rows, cols, ipiv, j, vn1, vn2);
                rankR++;
            }
            //return rankR;
        } 

        /// <summary>
        /// Find pivot column and check rank deficient
        /// </summary>
        /// <param name="ipiv"></param>
        /// <param name="vn1"></param>
        /// <param name="vn2"></param>
        /// <param name="pivotIndex"></param>
        /// <param name="criteria"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static bool PivotAndCheckRankDeficient(int[] ipiv, double[] vn1, double[] vn2, int pivotIndex, double criteria)
        {
            int n = ipiv.Length;

            int best = pivotIndex;
            double bestNorm = vn1[pivotIndex];

            for (int j = pivotIndex + 1; j < n; j++)
            {
                double v = vn1[j];
                if (v > bestNorm) { bestNorm = v; best = j; }
            }

            if (best != pivotIndex)
            {
                (ipiv[pivotIndex], ipiv[best]) = (ipiv[best], ipiv[pivotIndex]);
                (vn1[pivotIndex], vn1[best]) = (vn1[best], vn1[pivotIndex]);
                (vn2[pivotIndex], vn2[best]) = (vn2[best], vn2[pivotIndex]);
            }

            return bestNorm <= criteria;
        }

        /// <summary>
        /// Updating trailing column norm
        /// </summary>
        /// <param name="R"></param>
        /// <param name="ipiv"></param>
        /// <param name="pivot"></param>
        /// <param name="vn1"></param>
        /// <param name="vn2"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void UpdateTrailingColumnNorms(double[] arr, int rows, int cols, int[] ipiv, int pivot, double[] vn1, double[] vn2)
        {
 
            const double TOL3Z = 0.001;

            for (int pos = pivot + 1; pos < cols; pos++)
            {
                double v = vn1[pos];
                if (v <= EPS) continue;

                int col = ipiv[pos];
                int basec = col * rows;

                double r = Math.Abs(arr[basec + pivot]); // R[pivot, col] after update

                double ratio = r / v;
                double t = 1.0 - ratio * ratio;
                if (t < 0.0) t = 0.0;

                double newv = v * Math.Sqrt(t);
                vn1[pos] = newv;

                if (newv <= TOL3Z * vn2[pos])
                {
                    double s = 0.0;
                    for (int i = pivot + 1; i < rows; i++)
                    {
                        double e = arr[basec + i];
                        s += e * e;
                    }
                    double exact = Math.Sqrt(s);
                    vn1[pos] = exact;
                    vn2[pos] = exact;
                }
            }
        }


        /// <summary>
        /// Applies a Householder transformation to zero out subdiagonal elements of column pivot in A,
        /// and simultaneously applies the same transformation to vector b.
        /// Both A and b are modified in-place.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="ipiv"></param>
        /// <param name="pivot"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ApplyHouseholderToColumn(double[] arr, int rows, int cols, int[] ipiv, int pivot, Span<double> b)
        {
            int col = ipiv[pivot];
            int basev = col * rows;

            // norm2 of A[pivot:, col]
            double norm2 = 0.0;
            for (int i = pivot; i < rows; i++)
            {
                double e = arr[basev + i];
                norm2 += e * e;
            }
            if (norm2 <= EPS)
            {
                // keep triangular shape
                for (int i = pivot + 1; i < rows; i++) arr[basev + i] = 0.0;
                return;
            }

            double norm = Math.Sqrt(norm2);

            double a0 = arr[basev + pivot];
            double s = a0 >= 0.0 ? -1.0 : 1.0;
            double u1 = a0 - s * norm;
            if (Math.Abs(u1) <= EPS)
            {
                arr[basev + pivot] = s * norm;
                for (int i = pivot + 1; i < rows; i++) arr[basev + i] = 0.0;
                return;
            }

            double inv_u1 = 1.0 / u1;
            for (int i = pivot + 1; i < rows; i++) arr[basev + i] *= inv_u1;

            double tau = -s * u1 / norm;
            arr[basev + pivot] = s * norm;

            // Apply to remaining columns.
            // This is the dominant repeated kernel in QR.
            // The trailing panel shrinks as pivot increases, so switch dynamically:
            //   large trailing panel -> Parallel.For
            //   small trailing panel -> sequential for-loop
            ApplyHouseholderTrailingParallel(arr, rows, cols, ipiv, pivot, basev, tau);

            // Apply to b (sequential)
            double dotb = b[pivot];
            for (int i = pivot + 1; i < rows; i++)
                dotb += arr[basev + i] * b[i];

            double scale_b = tau * dotb;
            b[pivot] -= scale_b;
            for (int i = pivot + 1; i < rows; i++)
                b[i] -= scale_b * arr[basev + i];

            // store only R: zero out below diagonal in pivot column
            for (int i = pivot + 1; i < rows; i++)
                arr[basev + i] = 0.0;
        }

        private static void ApplyHouseholderTrailingParallel(
            double[] arr,
            int rows,
            int cols,
            int[] ipiv,
            int pivot,
            int basev,
            double tau)
        {
            Parallel.For(pivot + 1, cols, k =>
            {
                int basej = ipiv[k] * rows;

                double dot = arr[basej + pivot];
                for (int i = pivot + 1; i < rows; i++)
                    dot += arr[basev + i] * arr[basej + i];

                double scale = tau * dot;

                arr[basej + pivot] -= scale;
                for (int i = pivot + 1; i < rows; i++)
                    arr[basej + i] -= scale * arr[basev + i];
            });
        }
#endregion

#region Backward substitution i.e. inv(R11) Q^t b

        /// <summary>
        /// Backward substitution for the upper-triangular system R11 x1 = Q^T b.
        /// </summary>
        /// <param name="flatten_matrix"></param>
        /// <param name="rows"></param>
        /// <param name="Qtb"></param>
        /// <param name="x"></param>
        /// <param name="base_rank"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void BackwardSubstitutionForRHS(double[] flatten_matrix, int[] base_rank, double[] Qtb, double[] x)
        {          
            int r = base_rank.Length;
            Array.Copy(Qtb, x, r); // x1 = Q^T b (pivot order)

            for (int i = r - 1; i >= 0; i--)
            {
                double sum = x[i];

                for (int j = i + 1; j < r; j++)
                    sum -= flatten_matrix[base_rank[j] + i] * x[j];

                x[i] = sum / flatten_matrix[base_rank[i] + i];
            }
        }
#endregion

#region Under-determined / rank-deficient minimum 2-norm completion
        /// <summary>
        /// Completes the pivot-ordered solution to the minimum 2-norm solution
        /// for rank-deficient or under-determined cases.
        ///
        /// R is assumed to contain [R11 R12] in pivoted column order.
        /// The top rank rows of the free columns are overwritten by
        /// Y = R11^{-1} R12.
        ///
        /// x is stored in pivot order and modified in-place.
        /// </summary>
        private static bool CompleteMinimumNormSolution(
            double[] flatten_matrix,
            int rows,
            int cols,
            int[] ipiv,
            int[] base_rank,
            double[] x)
        {

            int r = base_rank.Length;
            int d = cols - r;

            if (r == 0)
            {
                // No range component; the minimum-norm solution is x = 0.
                Array.Clear(x, 0, cols);
                return true;
            }

            if (d <= 0)
            {
                // No free variables.
                return true;
            }

            ComputeYInPlace(flatten_matrix, rows, cols, r, ipiv, base_rank);

            // base_free[j] = base offset of pivot-free columns in original storage.
            int[] base_free = new int[d];
            for (int j = 0; j < d; j++)
                base_free[j] = ipiv[r + j] * rows;

            // Choose the smaller SPD system and solve it by dense Cholesky.
            if (d <= r)
            {
                return ApplyMinimumNormByYtYPlusICholesky(flatten_matrix, r, d, base_free, x);
            }
            else
            {
                return ApplyMinimumNormByYYtPlusICholesky(flatten_matrix, r, d, base_free, x);
            }
        }

        /// <summary>
        /// Computes Y = R11^{-1} R12 in-place.
        /// The top r rows of the free columns are overwritten by Y.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void ComputeYInPlace(
            double[] arr,
            int rows,
            int cols,
            int r,
            int[] pivot,
            int[] base_rank)
        {

            Parallel.For(r, cols, k =>
            {
                int colk = pivot[k];
                int basek = colk * rows;

                for (int i = r - 1; i >= 0; i--)
                {
                    double v = arr[basek + i];

                    for (int j = i + 1; j < r; j++)
                        v -= arr[base_rank[j] + i] * arr[basek + j];

                    v /= arr[base_rank[i] + i];
                    arr[basek + i] = v;
                }
            });

        }

        /// <summary>
        /// d <= r route.
        /// Solves (I + Y^T Y) x2 = Y^T x1 by dense Cholesky,
        /// then updates x1 <- x1 - Y x2.
        /// x is stored in pivot order and modified in-place.
        /// </summary>
        private static bool ApplyMinimumNormByYtYPlusICholesky(
            double[] arr,
            int r,
            int d,
            int[] base_free,
            double[] x)
        {
            double[] spd = new double[d * (d + 1) / 2]; // packed lower storage
            double[] rhs = new double[d];

            BuildYtYPlusI(arr, r, d, base_free, spd);

            // rhs = Y^T x1
            for (int j = 0; j < d; j++)
            {
                double sum = 0.0;
                int basec = base_free[j];

                for (int i = 0; i < r; i++)
                    sum += arr[basec + i] * x[i];

                rhs[j] = sum;
            }

            if (!CholeskyDecomposeLower(spd, d)) return false;
            CholeskySolveLowerInPlace(spd, rhs, d); // rhs becomes x2

            // x1 <- x1 - Y x2
            for (int i = 0; i < r; i++)
            {
                double dot = 0.0;

                for (int j = 0; j < d; j++)
                    dot += arr[base_free[j] + i] * rhs[j];

                x[i] -= dot;
            }

            // x2 write-back.
            for (int j = 0; j < d; j++)
                x[r + j] = rhs[j];

            return true;
        }

        /// <summary>
        /// d > r route.
        /// Solves (I + Y Y^T) z = x1 by dense Cholesky,
        /// then sets x1 <- z and x2 <- Y^T z.
        /// x is stored in pivot order and modified in-place.
        /// </summary>
        private static bool ApplyMinimumNormByYYtPlusICholesky(
            double[] arr,
            int r,
            int d,
            int[] base_free,
            double[] x)
        {
            double[] spd = new double[r * (r + 1) / 2]; // packed lower storage
            double[] z = new double[r];

            BuildYYtPlusI(arr, r, d, base_free, spd);

            Array.Copy(x, z, r);

            if (!CholeskyDecomposeLower(spd, r)) return false;
            CholeskySolveLowerInPlace(spd, z, r); // z becomes (I + Y Y^T)^-1 x1

            // x1 <- z
            for (int i = 0; i < r; i++)
                x[i] = z[i];

            // x2 <- Y^T z
            for (int j = 0; j < d; j++)
            {
                double sum = 0.0;
                int basec = base_free[j];

                for (int i = 0; i < r; i++)
                    sum += arr[basec + i] * z[i];

                x[r + j] = sum;
            }

            return true;
        }

        /// <summary>
        /// Packed lower-triangular index for row >= col.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int PackedLowerIndex(int row, int col)
        {
            return row * (row + 1) / 2 + col;
        }

        /// <summary>
        /// Builds S = I + Y^T Y in packed lower-triangular storage.
        /// Entry S[row, col] with row >= col is stored at row * (row + 1) / 2 + col.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void BuildYtYPlusI(double[] arr, int r, int d, int[] base_free, double[] spd)
        {
            //for (int row = 0; row < d; row++)
            Parallel.For(0, d, row =>
            {
                int baserow = base_free[row];

                for (int col = 0; col <= row; col++)
                {
                    int basecol = base_free[col];
                    double sum = (row == col) ? 1.0 : 0.0;

                    for (int i = 0; i < r; i++)
                        sum += arr[baserow + i] * arr[basecol + i];

                    spd[PackedLowerIndex(row, col)] = sum;
                }
            }
            );
        }

        /// <summary>
        /// Builds S = I + Y Y^T in packed lower-triangular storage.
        /// Entry S[row, col] with row >= col is stored at row * (row + 1) / 2 + col.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void BuildYYtPlusI(double[] arr, int r, int d, int[] base_free, double[] spd)
        {
            //for (int row = 0; row < r; row++)
            Parallel.For(0, r, row =>
            {
                for (int col = 0; col <= row; col++)
                {
                    double sum = (row == col) ? 1.0 : 0.0;

                    for (int j = 0; j < d; j++)
                        sum += arr[base_free[j] + row] * arr[base_free[j] + col];

                    spd[PackedLowerIndex(row, col)] = sum;
                }
            }
            );
        }

        /// <summary>
        /// In-place Cholesky factorization A = L L^T for packed lower storage.
        /// The packed array is overwritten by the packed lower Cholesky factor L.
        /// </summary>
        private static bool CholeskyDecomposeLower(double[] a, int n)
        {
            for (int i = 0; i < n; i++)
            {
                int rowiBase = i * (i + 1) / 2;

                for (int j = 0; j <= i; j++)
                {
                    double sum = a[rowiBase + j];
                    int rowjBase = j * (j + 1) / 2;

                    for (int k = 0; k < j; k++)
                        sum -= a[rowiBase + k] * a[rowjBase + k];

                    if (i == j)
                    {
                        if (sum <= 0.0 || double.IsNaN(sum))
                            return false;

                        a[rowiBase + j] = Math.Sqrt(sum);
                    }
                    else
                    {
                        a[rowiBase + j] = sum / a[rowjBase + j];
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Solves L L^T x = b in-place for packed lower Cholesky factor L.
        /// The vector b is overwritten by x.
        /// </summary>
        private static void CholeskySolveLowerInPlace(double[] lower, double[] b, int n)
        {
            // Forward solve: L y = b.
            for (int i = 0; i < n; i++)
            {
                double sum = b[i];
                int rowiBase = i * (i + 1) / 2;

                for (int j = 0; j < i; j++)
                    sum -= lower[rowiBase + j] * b[j];

                b[i] = sum / lower[rowiBase + i];
            }

            // Backward solve: L^T x = y.
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = b[i];

                for (int j = i + 1; j < n; j++)
                    sum -= lower[PackedLowerIndex(j, i)] * b[j];

                b[i] = sum / lower[PackedLowerIndex(i, i)];
            }
        }

        /// <summary>
        /// Compute residual (||Ax - b|| = ||RPx - Qtb||)
        /// </summary>
        /// <param name="R"></param>
        /// <param name="Qtb"></param>
        /// <param name="ipiv"></param>
        /// <param name="perm_x"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double GetResidualNorm(ReadOnlySpan<double> arr, int rows, int cols, ReadOnlySpan<double> Qtb, int[] ipiv, ReadOnlySpan<double> perm_x) // unpivot前の x
        {

            int k = Math.Min(rows, cols);

            double s2 = 0.0;

            // head: r_i = Qtb[i] - sum_{j=i..n-1} R[i, ipiv[j]] * x[j]
            for (int i = 0; i < k; i++)
            {
                double ri = Qtb[i];
                for (int j = i; j < cols; j++)
                    ri -= arr[ipiv[j] * rows + i] * perm_x[j];

                s2 += ri * ri;
            }

            // tail: Qtb[k:rows]
            for (int i = k; i < rows; i++)
            {
                double e = Qtb[i];
                s2 += e * e;
            }

            return Math.Sqrt(s2);
        }
#endregion

#region  Unpivot solution
        /// <summary>
        /// Reverses the column pivoting applied to the solution vector.
        /// </summary>
        /// <param name="cols"></param>
        /// <param name="ipiv"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        static double[] UnpivotSolution(int cols, int[] ipiv, double[] x)
        {
            double[] solution = new double[cols];
            for (int j = 0; j < cols; j++) solution[ipiv[j]] = x[j];
            return solution;
        }
#endregion

#region Public result object and status codes
    public enum LSQSolverStatus
    {
        Success = 0,
        NullMatrix,
        EmptyMatrix,
        NullVector,
        DimensionMismatch,
        CholeskyFailed,
            InvalidVector,
            InvalidMatrix
        }

    /// <summary>
    /// This object contains the result of LSQSolver.Solve().
    /// </summary>
    public class LSQSolverResult
    {

        public object? Tag { get; internal set; } = null;

        /// <summary>
        /// Status code of the solve operation.
        /// </summary>
        public LSQSolverStatus Status { get; internal set; } = LSQSolverStatus.Success;
 
        /// <summary>
        /// Solution of LSQ problem.
        /// </summary>
        public double[] Solution { get; internal set; } = Array.Empty<double>();

        /// <summary>
        /// Rows of A.
        /// </summary>
        public int Rows { get; internal set; } = 0;

        /// <summary>
        /// Cols of A.
        /// </summary>
        public int Cols { get; internal set; } = 0;

        /// <summary>
        /// Estimated rank of R.
        /// </summary>
        public int Rank { get; internal set; } = 0;

        /// <summary>
        /// Residual norm.
        /// </summary>
        public double ResidualNorm { get; internal set; } = 0;

        /// <summary>
        /// R factor after column-pivoted QR.
        /// Columns are stored in pivoted order.
        /// </summary>
        public MatrixObject? R { get; internal set; } = null;

        /// <summary>
        /// Transformed right-hand side Q^T b.
        /// </summary>
        public double[]? Qtb { get; internal set; } = null;

        /// <summary>
        /// Pivot[j] gives the original column index of the j-th pivoted column.
        /// </summary>
        public int[]? Pivot { get; internal set; } = null;

    

        /// <summary>
        /// Convert current state to string.
        /// </summary>
        /// <returns>state</returns>
        public string ToString(bool omit = false)
        {
            string ret = string.Empty;
            int displayMax;

            ret += "R:\n";
            if (R is null)
            {
                ret += "<not stored>\n";
            }
            else
            {
                for (int i = 0; i < R.Rows; i++)
                {
                    if (omit && i >= 10) break;
                    for (int j = 0; j < R.Cols; j++)
                    {
                        if (omit && j >= 10) break;
                        ret += $"{R[i, j]:F6}\t";
                    }
                    if (omit && R.Cols > 10) ret += "...";
                    ret += "\n";
                }
                if (omit && R.Rows > 10) ret += "...\n";
            }
            ret += "\n";

            ret += "Qtb:\n";
            if (Qtb is null)
            {
                ret += "<not stored>\n";
            }
            else if (omit)
            {
                displayMax = Math.Min(Qtb.Length, 10);
                for (int i = 0; i < displayMax; i++)
                {
                    if (i > 0) ret += "\t";
                    ret += $"{Qtb[i]:F6}";
                }

                if (Qtb.Length > displayMax) ret += "...";
                ret += "\n";
            }
            else
            {
                for (int i = 0; i < Qtb.Length; i++)
                {
                    if (i > 0) ret += "\t";
                    ret += $"{Qtb[i]:F6}";
                }
                ret += "\n";
            }
            ret += "\n";

            ret += "Pivot:\n";
            if (Pivot is null)
            {
                ret += "<not stored>\n";
            }
            else if (omit)
            {
                displayMax = int.Min(Pivot.Length, 10);
                ret += string.Join("\t", Pivot[0..displayMax]);
                if (Pivot.Length > displayMax) ret += "...";
                ret += "\n";
            }
            else
            {
                ret += string.Join(" ", Pivot) + "\n";
            }
            ret += "\n";

            ret += $"Status:\n{Status}\n\n";
            ret += $"Rows:\n{Rows}\n\n";
            ret += $"Cols:\n{Cols}\n\n";
            ret += $"Rank:\n{Rank}\n\n";

            ret += "x:\n";
            if (omit)
            {
                displayMax = Math.Min(Solution.Length, 10);
                for (int i = 0; i < displayMax; i++)
                {
                    if (i > 0) ret += "\t";
                    ret += $"{Solution[i]:F6}";
                }

                if (Solution.Length > displayMax) ret += "...";
                ret += "\n";
            }
            else
            {
                for (int i = 0; i < Solution.Length; i++)
                {
                    if (i > 0) ret += "\t";
                    ret += $"{Solution[i]:F6}";
                }
                ret += "\n";
            }
            ret += "\n";

            ret += $"Residual norm:\n{ResidualNorm}\n\n";

            return ret;
        }
    }

    }
#endregion

#region MatrixObject

    /// <summary>
    /// Matrix object for LSQ solver, storing values in column-major order.
    /// </summary>
    /// <remarks>It stores double-precision matrix data as COLUMN-major 1D array.</remarks>
    public class MatrixObject
    {
        /// <summary>
        /// internal storage
        /// </summary>
        internal readonly double[] array;
        /// <summary>
        /// Number of rows
        /// </summary>
        public int Rows { get; }
        /// <summary>
        /// Number of columns
        /// </summary>
        public int Cols { get; }
        /// <summary>
        /// Constructor with memory allocation
        /// </summary>
        /// <param name="rows"></param>
        /// <param name="cols"></param>
        public MatrixObject(int rows, int cols)
        {
            Rows = rows;
            Cols = cols;
            array = new double[rows * cols];
        }
        /// <summary>
        /// Copy constructor
        /// </summary>
        /// <param name="other"></param>
        public MatrixObject(MatrixObject other)
        {
            Rows = other.Rows;
            Cols = other.Cols;
            array = (double[])other.array.Clone();
        }
        /// <summary>
        /// Copy constructor with pivot
        /// </summary>
        /// <param name="other"></param>
        public MatrixObject(MatrixObject other, int[] pivot)
        {
            Rows = other.Rows;
            Cols = other.Cols;

            array = new double[other.array.Length];
            for (int j = 0; j < Cols; j++)
            {
                for (int i = 0; i < Rows; i++)
                {
                    array[i + Rows * j] = other[i, pivot[j]];
                }
            }
        }
        /// <summary>
        /// Constructor by jagged array
        /// </summary>
        /// <param name="A"></param>
        /// <exception cref="ArgumentException"></exception>
        public MatrixObject(double[][] A)
        {
            if (A == null || A.Length == 0 || A[0] == null)
                throw new ArgumentException("Jagged array must be non-null and non-empty.");

            Rows = A.Length;
            Cols = A[0].Length;
            array = new double[Rows * Cols];

            for (int i = 0; i < Rows; i++)
            {
                if (A[i].Length != Cols)
                    throw new ArgumentException("All rows in the A array must have the same length.");

                for (int j = 0; j < Cols; j++)
                    array[i + Rows * j] = A[i][j];
            }
        }
        /// <summary>
        /// Index access for element
        /// </summary>
        /// <param name="i">i-th row</param>
        /// <param name="j">j-th column</param>
        /// <returns></returns>
        public ref double this[int i, int j]
        {
            get => ref array[i + Rows * j];
        }
    }

#endregion

}

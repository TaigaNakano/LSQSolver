# Polynomial Fitting

## Overview

This page explains the basic usage of `LSQSolver` using polynomial fitting as an example. Polynomial fitting is the problem of finding a line or curve that passes near a given set of data points. For example, suppose we have the following data.

```text
x: 0, 1, 2, 3
y: 1.0, 2.1, 2.9, 4.2
```

We can fit a line, a quadratic curve, or another polynomial to these data points. The scatter plot of the data is shown below.

![](./figs/polynomial_fit/original_scatter.svg)

Note that `LSQSolver` is not a library specialized for polynomial fitting. It is a library for solving linear least-squares problems represented by a matrix and a vector. In this page, polynomial fitting is used as an example of a problem that can be rewritten as a linear least-squares problem handled by `LSQSolver`.

---

## Formulating polynomial fitting

Let us consider fitting a polynomial to data points. We write the data points as $`(x_1, y_1), (x_2, y_2), \ldots, (x_m, y_m)`$, where $`m`$ is the number of data points, $`x_i`$ is the horizontal coordinate of the $`i`$-th point, and $`y_i`$ is the vertical coordinate of the $`i`$-th point.

For example, if we use a quadratic polynomial, it has the following form.

```math
p(x) = c_0 + c_1 x + c_2 x^2.
```

Here, $`p(x)`$ is the value of the curve at $`x`$, and $`c_0, c_1, c_2`$ are the coefficients that determine the shape of the polynomial. These coefficients are the unknowns we want to compute. Changing the coefficients changes the shape of the curve. In other words, polynomial fitting means choosing the shape of the curve so that it fits the given data points as well as possible.

If the points are roughly aligned along a straight line, a linear polynomial may be enough. If the points follow a slightly curved trend, a quadratic or cubic polynomial may fit them more naturally.

![](./figs/polynomial_fit/degree1_vs_degree2.svg)

The important point is that `LSQSolver` does not treat polynomials as a special object. If we regard the coefficients $`c_0, c_1, c_2`$ as unknowns, the problem can be written as a linear problem using a matrix and a vector. `LSQSolver` solves that linear problem and returns the coefficients.

In C#, we first build the matrix `A` corresponding to the polynomial coefficients, and then pass it to `LSQSolver.Solve()`. For a quadratic polynomial, the computation can be written as follows.

```csharp
double[] xs = { 0.0, 1.0, 2.0, 3.0 };
double[] ys = { 1.0, 2.1, 2.9, 4.2 };

// p(x) = c0 + c1 x + c2 x^2
double[][] data =
{
    new double[] { 1.0, xs[0], xs[0] * xs[0] },
    new double[] { 1.0, xs[1], xs[1] * xs[1] },
    new double[] { 1.0, xs[2], xs[2] * xs[2] },
    new double[] { 1.0, xs[3], xs[3] * xs[3] }
};

var A = new MatrixObject(data);
var result = Solve(A, ys, overwrite: false);

double[] c = result.Solution;

Console.WriteLine($"p(x) = {c[0]} + {c[1]} x + {c[2]} x^2");
Console.WriteLine($"residual = {result.ResidualNorm}");
```

Output:

```console
p(x) = 1.0399999999999996 + 0.890000000000001 x + 0.04999999999999968 x^2
residual = 0.17888543819998348
```

In this code, `xs` stores the horizontal coordinates and `ys` stores the vertical coordinates. The array `result.Solution` contains the polynomial coefficients $`c_0, c_1, c_2`$.

Depending on the number of points, the degree of the polynomial, and possible duplication in the data, the following situations can occur.

1. The polynomial is uniquely determined.
2. No polynomial of the chosen degree can pass through all points exactly.
3. Many polynomials satisfy the given conditions.
4. The rank drops.

In the following sections, we first explain situations 1–3 using figures and code, without assuming much mathematical background. After that, we summarize the same ideas using linear algebra and then discuss situation 4, where the rank drops.

---

## Situation 1: The polynomial is uniquely determined

First, consider the case where the number of data points matches the number of unknown coefficients. For example, the quadratic polynomial

```math
p(x) = c_0 + c_1 x + c_2 x^2
```

has three unknown coefficients. If we have three data points, and the points are in a suitable configuration, we can determine a quadratic polynomial that passes through all three points. As an example, consider the following three points.

```math
(-1, 1),\quad (0, 0),\quad (1, 1).
```

The quadratic polynomial passing through these points is

```math
p(x) = x^2.
```

Indeed, the following equalities hold.

```math
p(-1) = 1,\qquad p(0) = 0,\qquad p(1) = 1.
```

![](./figs/polynomial_fit/well_defined.svg)

The following code computes this result using `LSQSolver`.

```csharp
double[] xs = { -1.0, 0.0, 1.0 };
double[] ys = {  1.0, 0.0, 1.0 };

// p(x) = c0 + c1 x + c2 x^2
double[][] data =
{
    new double[] { 1.0, xs[0], xs[0] * xs[0] },
    new double[] { 1.0, xs[1], xs[1] * xs[1] },
    new double[] { 1.0, xs[2], xs[2] * xs[2] }
};

var A = new MatrixObject(data);
var result = Solve(A, ys, overwrite: false);

double[] c = result.Solution;

Console.WriteLine($"p(x) = {c[0]} + {c[1]} x + {c[2]} x^2");
Console.WriteLine($"residual = {result.ResidualNorm}");
```

Output:

```console
p(x) = -0 + -0 x + 1 x^2
residual = 0
```

In this case, `residual` is zero. This means that the computed curve passes through the given points exactly.

---

## Situation 2: There are too many points to pass through exactly

Next, consider the case where there are many data points. For example, suppose we try to pass a line

```math
p(x) = c_0 + c_1 x
```

through the following four points.

```math
(0, 1.0),\quad (1, 2.1),\quad (2, 2.9),\quad (3, 4.2).
```

These four points are roughly aligned along a line. However, they are not exactly on a single straight line. Therefore, there is no line that passes through all four points exactly. Even in this case, we can still find a line that fits the data well overall.

![](./figs/polynomial_fit/over_determined.svg)

The following code computes such a line using `LSQSolver`.

```csharp
double[] xs = { 0.0, 1.0, 2.0, 3.0 };
double[] ys = { 1.0, 2.1, 2.9, 4.2 };

// p(x) = c0 + c1 x
double[][] data =
{
    new double[] { 1.0, xs[0] },
    new double[] { 1.0, xs[1] },
    new double[] { 1.0, xs[2] },
    new double[] { 1.0, xs[3] }
};

var A = new MatrixObject(data);
var result = Solve(A, ys, overwrite: false);

double[] c = result.Solution;

Console.WriteLine($"p(x) = {c[0]} + {c[1]} x");
Console.WriteLine($"residual = {result.ResidualNorm}");
```

Output:

```console
p(x) = 0.9900000000000004 + 1.0399999999999998 x
residual = 0.2049390153191921
```

In this case, `residual` is not zero. This means that the line does not pass through all points exactly. Even so, `LSQSolver` returns a line that is close to the data points overall.

---

## Situation 3: Many polynomials satisfy the data

Now consider the case where there are too many unknowns compared with the number of data points. Suppose we are given only two points.

```math
(0, 1),\quad (1, 2).
```

There is only one line passing through these two points. However, if we try to use a cubic polynomial

```math
p(x) = c_0 + c_1 x + c_2 x^2 + c_3 x^3
```

then the answer is no longer unique. There are many cubic curves passing through the same two points.

![](./figs/polynomial_fit/under_determined.svg)

The following code computes one such cubic polynomial using `LSQSolver`.

```csharp
double[] xs = { 0.0, 1.0 };
double[] ys = { 1.0, 2.0 };

// p(x) = c0 + c1 x + c2 x^2 + c3 x^3
double[][] data =
{
    new double[] { 1.0, xs[0], xs[0] * xs[0], xs[0] * xs[0] * xs[0] },
    new double[] { 1.0, xs[1], xs[1] * xs[1], xs[1] * xs[1] * xs[1] }
};

var A = new MatrixObject(data);
var result = Solve(A, ys, overwrite: false);

double[] c = result.Solution;

Console.WriteLine($"p(x) = {c[0]} + {c[1]} x + {c[2]} x^2 + {c[3]} x^3");
Console.WriteLine($"residual = {result.ResidualNorm}");
```

Output:

```console
p(x) = 0.9999999999999999 + 0.33333333333333326 x + 0.3333333333333332 x^2 + 0.3333333333333335 x^3
residual = 0
```

In this case, there are many curves passing through the data points. `LSQSolver` chooses a natural answer among them, avoiding unnecessarily large coefficients. The important point is that when multiple answers explain the data, `LSQSolver` returns one manageable solution.

---

## What `LSQSolver` returns

So far, we have seen three representative situations that can occur in polynomial fitting. If the number of points and the number of coefficients match appropriately, a curve passing exactly through the points may be obtained. If there are too many points, the curve may not pass through all points exactly, so we look for a curve that fits well overall. If there are too many coefficients, many curves may satisfy the data, so one natural curve is selected.

The following figures show the polynomials obtained using `LSQSolver` for the three cases.

<p align="center">
  <img src="./figs/polynomial_fit/well_defined.svg" width="30%" />
  <img src="./figs/polynomial_fit/over_determined.svg" width="30%" />
  <img src="./figs/polynomial_fit/under_determined.svg" width="30%" />
</p>

At this point, the basic idea of using `LSQSolver` for polynomial fitting has been explained. In practice, we build the matrix `A` from the data points, call `Solve()`, and read the coefficients from the result. On the other hand, what do phrases such as “fits well overall” and “a natural curve” mean in linear algebra? The following sections explain these meanings in more detail.

---

## Polynomial fitting as a system of linear equations

From here, we describe why polynomial fitting can be solved by `LSQSolver` using linear algebra. The previous sections already explain the basic usage. The following sections are supplementary material for readers who want to understand more precisely what kind of solution `LSQSolver` returns.

Let the data points be $`(x_1, y_1), (x_2, y_2), \ldots, (x_m, y_m)`$, and let the degree-$`n`$ polynomial be

```math
p(x) = c_0 + c_1 x + c_2 x^2 + \cdots + c_n x^n.
```

Here, $`m`$ is the number of data points, $`n`$ is the degree of the polynomial, and $`c_0, c_1, \ldots, c_n`$ are the coefficients to be computed. The number of coefficients is $`n+1`$. For each data point, we want the following relation to hold approximately.

```math
p(x_i) \approx y_i
\qquad (i = 1,2,\ldots,m)
```

This can be written in matrix form as

```math
A\mathbf{c} \approx \mathbf{y}.
```

Here, $`A`$, $`\mathbf{c}`$, and $`\mathbf{y}`$ are defined as follows.

```math
A =
\begin{bmatrix}
1 & x_1 & x_1^2 & \cdots & x_1^n \\
1 & x_2 & x_2^2 & \cdots & x_2^n \\
\vdots & \vdots & \vdots & & \vdots \\
1 & x_m & x_m^2 & \cdots & x_m^n
\end{bmatrix},
\qquad
\mathbf{c} =
\begin{bmatrix}
c_0 \\
c_1 \\
\vdots \\
c_n
\end{bmatrix},
\qquad
\mathbf{y} =
\begin{bmatrix}
y_1 \\
y_2 \\
\vdots \\
y_m
\end{bmatrix}.
```

The matrix $`A`$ is an $`m \times (n+1)`$ matrix. The number of rows $`m`$ corresponds to the number of data points, and the number of columns $`n+1`$ corresponds to the number of polynomial coefficients. In this way, polynomial fitting becomes a linear algebra problem for computing the coefficient vector $`\mathbf{c}`$.

---

## Matrix shape and rank

The nature of the problem can be organized by the shape and rank of the matrix $`A`$. Let $`r = \mbox{rank}(A)`$ be the rank of $`A`$. Here, $`r`$ represents the number of independent pieces of information contained in the matrix $`A`$.

| Situation             |   Example condition | Nature of the solution                  |
| --------------------- | ------------------: | --------------------------------------- |
| Exactly determined    | $`m = n+1,\ r = n+1`$ | Unique solution                         |
| Too many points       |           $`m > n+1`$ | In general, cannot be satisfied exactly |
| Too many coefficients |           $`m < n+1`$ | In general, the solution is not unique  |
| Rank deficient        |   $`r < \min(m,n+1)`$ | Independent information is missing      |

`LSQSolver` estimates the numerical rank. Therefore, it considers not only the apparent matrix size, but also how much independent information is actually contained in the matrix.

---

## Least-squares solution

Consider the case where the number of data points is larger than the number of coefficients, that is, $`m > n+1`$. In this case, the equation $`A\mathbf{c} = \mathbf{y}`$ generally cannot be satisfied exactly. In terms of polynomial fitting, this means that there is no curve of the chosen degree that passes through all data points exactly. Instead of requiring exact agreement, we choose the coefficient vector $`\mathbf{c}`$ that makes the residual $`A\mathbf{c} - \mathbf{y}`$ as small as possible.

Such a solution is called a **least-squares solution**. It is the solution of the following minimization problem.

```math
\min_{\mathbf{c}} \|A\mathbf{c} - \mathbf{y}\|_2.
```

Here, $`A\mathbf{c} - \mathbf{y}`$ is the residual vector, and $`\|A\mathbf{c} - \mathbf{y}\|_2`$ is the size of the residual. The notation $`|\cdot|_2`$ denotes the Euclidean norm, which is the square root of the sum of squares of the components. In the context of polynomial fitting, this means choosing a curve that fits the data well overall, even when it cannot pass through all points exactly.

The phrase “a curve that fits well overall” used earlier corresponds to this least-squares solution. In the earlier example with too many points, a single line cannot pass through all four data points exactly. Therefore, `LSQSolver` returns coefficients that make the residual small.

---

## Minimum-norm solution

On the other hand, consider the case where the number of coefficients is larger than the number of conditions, that is, $`m < n+1`$. In this case, there may be multiple vectors $`\mathbf{c}`$ satisfying $`A\mathbf{c} = \mathbf{y}`$. In other words, multiple polynomials may explain the same data points. We then need a criterion for choosing one solution.

In such a case, `LSQSolver` returns a **minimum-norm solution**. A minimum-norm solution is the solution whose coefficient vector has the smallest size $`\|\mathbf{c}\|_2`$ among all solutions satisfying the conditions. It can be written as the following constrained minimization problem.

```math
\min_{\mathbf{c}} \|\mathbf{c}\|_2
\quad
\text{subject to}
\quad
A\mathbf{c} = \mathbf{y}.
```

The phrase “a natural curve” used earlier corresponds to this minimum-norm solution. Here, “natural” is not determined only by the visual appearance of the curve. It is defined by the linear-algebraic criterion that the coefficient vector has the smallest norm $`\|\mathbf{c}\|_2`$.

For example, there are infinitely many cubic curves passing through two points. Among them, `LSQSolver` chooses a solution that satisfies the given conditions while avoiding unnecessarily large coefficients. In this sense, the earlier phrase “a manageable solution” corresponds to the minimum-norm solution.

However, “a natural curve” does not necessarily mean the visually smoothest curve. In this document, it means the solution whose coefficient vector has the smallest norm.

---

## Situation 4: The rank drops

In polynomial fitting, the rank may drop. For example, this can happen when the same $`x`$ value appears more than once.

```text
x: 0, 0, 1
y: 1, 1, 2
```

Even if we try to fit a quadratic polynomial, the first and second data points represent the same condition. In this case, the matrix $`A`$ becomes

```math
A =
\begin{bmatrix}
1 & 0 & 0^2 \\
1 & 0 & 0^2 \\
1 & 1 & 1^2
\end{bmatrix}
=
\begin{bmatrix}
1 & 0 & 0 \\
1 & 0 & 0 \\
1 & 1 & 1
\end{bmatrix}.
```

The first and second rows are identical, so there are only two linearly independent row vectors. In this document, we say that the rank drops when the number of independent conditions is smaller than the apparent number of conditions.

Next, consider the case where slightly different $`y`$ values are given for the same $`x`$ value.

```text
x: 0, 0, 1
y: 1, 1.25, 2
```

In this case, two values, $`y=1`$ and $`y=1.25`$, are given at $`x=0`$. Therefore, no single polynomial can satisfy both conditions exactly. Let us compute this case using `LSQSolver`.

```csharp
double[] xs = { 0.0, 0.0, 1.0 };
double[] ys = { 1.0, 1.25, 2.0 };

// p(x) = c0 + c1 x + c2 x^2
double[][] data =
{
    new double[] { 1.0, xs[0], xs[0] * xs[0] },
    new double[] { 1.0, xs[1], xs[1] * xs[1] },
    new double[] { 1.0, xs[2], xs[2] * xs[2] }
};

var A = new MatrixObject(data);
var result = Solve(A, ys, overwrite: false);

double[] c = result.Solution;

Console.WriteLine($"p(x) = {c[0]} + {c[1]} x + {c[2]} x^2");
Console.WriteLine($"rank = {result.Rank}");
Console.WriteLine($"residual = {result.ResidualNorm}");
```

Output:

```console
p(x) = 1.1249999999999998 + 0.43750000000000017 x + 0.43750000000000006 x^2
rank = 2
residual = 0.1767766952966373
```

![](./figs/polynomial_fit/rank_deficient.svg)

In this case, the matrix size appears to indicate three conditions. However, the number of independent conditions is actually smaller. `LSQSolver` estimates the numerical rank and returns a solution according to that rank, even for such rank-deficient problems.

When the rank drops, the three situations described earlier may no longer be cleanly separated. Since independent information is missing, the solution may fail to be unique. In such a case, `LSQSolver` returns a minimum-norm solution, choosing the candidate whose coefficient vector has the smallest norm.

On the other hand, if the conditions are inconsistent, as in the case where different $`y`$ values are given for the same $`x`$, no curve can pass through all points exactly. In that case, `LSQSolver` first finds a curve that makes the residual small, that is, a least-squares solution that fits the data well overall. If multiple such solutions exist, it then chooses the minimum-norm solution among them.

---

## Summary

Polynomial fitting looks like the problem of fitting a curve to points. However, `LSQSolver` is not a library specialized for polynomial fitting. By treating the polynomial coefficients as unknowns, polynomial fitting can be written as the following linear problem and solved by `LSQSolver`.

```math
A\mathbf{c} \approx \mathbf{y}
```

Here, $`A`$ is the matrix built from the $`x`$ values of the data points, $`\mathbf{c}`$ is the vector of polynomial coefficients, and $`\mathbf{y}`$ is the vector of the $`y`$ values.

Depending on the number of points, the degree of the polynomial, and duplication in the data, the problem may fall into one of the following situations. The polynomial passing through the given points may be

* uniquely determined,
* unable to pass through all points exactly,
* non-unique, or
* rank deficient.

In linear algebra terms, these situations can be organized using rank, least-squares solutions, and minimum-norm solutions.

* A curve that “fits well overall” corresponds to a **least-squares solution**, which makes the residual small.
* A “natural curve” corresponds to a **minimum-norm solution**, which makes the coefficient vector small.
* When the rank drops, independent information is missing, and both least-squares and minimum-norm ideas may be involved.

In this way, polynomial fitting is a natural example for explaining the basic usage and features of `LSQSolver`.

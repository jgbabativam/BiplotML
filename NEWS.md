# BiplotML 1.1.1

## Resubmission to CRAN

This version restores the package to CRAN after it was archived on 2023-10-29
due to a dependency on the `optimr` package, which was removed from CRAN at the
maintainer's request. All calls to `optimr()` now use `optimx::optimr()`, which
is the current home of that function.

### Dependency changes

* The `optimr` package has been removed from `Imports`. Its function `optimr()`
  is now imported from `optimx`, which absorbed it.
* `optimr`, `optimx`, and `shapes` moved from `Depends` to `Imports`, in line
  with CRAN policy.
* `RSpectra` moved from `Suggests` to `Imports`, since it is unconditionally
  required for large matrices in `proj_LogBip()`.

### Bug fixes

* `bootBLB()`: fixed a column-naming bug that caused `plotBLB()` to fail with
  "Column 'Dim1' not found" after bootstrap aggregation.
* `bootBLB()`: suppressed verbose control messages printed to the console by
  `optimx::optimr()` when using the CG method.
* Fixed typo `Sensitivy` -> `Sensitivity` in the confusion-matrix output of
  `bootBLB()` and `pred_LB()`.

### Internal changes

* All non-ASCII characters (typographic dashes, accented letters) in R source
  files replaced with portable ASCII equivalents, as required by CRAN.
* `utils::globalVariables()` declaration moved to `R/zzz.R` to prevent a
  roxygen2 block error in `plotBLB.R`.
* Removed incorrect `importFrom(stats, svd)` and `importFrom(stats, eigen)`
  from NAMESPACE: both functions belong to base R.
* All `ggplot2::aes()` calls updated to use `.data[[]]` tidy evaluation,
  eliminating R CMD check NOTEs about undefined global variables.
* `inst/CITATION` updated from the deprecated `citEntry()` to `bibentry()`.
* `sdv_MM()`: loop variable `j` initialised before the loop to prevent a
  potential "object not found" error on zero iterations.
* `RoxygenNote` updated to 7.3.3.

---

# BiplotML 1.1.0

## New algorithm: projection-based logistic biplot with missing data

Version 1.1.0 introduced a major new fitting method for the logistic biplot
model, described in:

> Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2022).
> A coordinate descent MM algorithm for logistic biplot model with missing data.
> *In process*.

### The problem addressed

All previous logistic biplot algorithms (alternating, external logistic, and
the conjugate gradient / iterated SVD methods introduced in v1.0.0) share a
structural limitation: each row of the data matrix has its own parameter vector
**theta**_i = mu_i + sum_s a_is * b_s. Consequently:

* The number of parameters to estimate **grows with the number of rows**,
  making the model computationally impractical for large datasets.
* Projecting new individuals as supplementary rows requires **re-running the
  full optimisation**, because the row markers are free parameters. This
  introduces a risk of overfitting.
* None of the previous algorithms handled **missing data** directly.

### The new approach

The new method reformulates the logistic biplot using Pearson's (1901) data
projection idea, extended to the logistic case by Landgraf & Lee (2020). Instead
of treating each row's coordinates as independent free parameters, the row
markers are expressed as a projection of the (centred) data matrix onto a
low-rank subspace **V**:

    A = (X - 1 * mu') * V

This single change has three important consequences:

1. **The number of parameters no longer depends on n.** Only the p x k matrix
   **V** (column markers) and the p-vector **mu** (intercepts) need to be
   estimated, regardless of how many rows the data matrix has.

2. **New individuals can be projected without refitting.** Given estimated **V**
   and **mu**, the row markers of any new observation x_new are simply:
   a_new = (x_new - mu') * V. No optimisation is required.

3. **Missing data are handled natively.** A weight matrix **W** (W_ij = 1 if
   x_ij is observed, 0 if missing) is incorporated into the loss function.
   Missing entries are imputed at each iteration using the current fitted values
   and a per-variable threshold that minimises the Balanced Accuracy (BACC),
   ensuring that classification performance is optimised throughout fitting.

### The algorithm

The objective function -- the negative log-likelihood weighted by **W** -- is
non-convex. To avoid dealing with it directly, it is majorized at each iteration
by a quadratic surrogate (the MM step), following the approach of
Babativa-Marquez & Vicente-Villardon (2021). The surrogate is then minimised
using a block coordinate descent algorithm:

* **Update V:** Fix mu; the optimal V consists of the k leading eigenvectors of
  a symmetric matrix derived from the surrogate, computed via eigendecomposition
  (using `RSpectra::eigs_sym()` for large matrices).
* **Update mu:** Fix V; the optimal mu has a closed-form solution as a weighted
  column mean.
* **Impute missing values:** Replace missing entries with binary predictions
  from the current fitted model, using per-variable BACC-optimal thresholds.
* **Repeat** until the relative decrease in the loss function falls below the
  convergence tolerance (default 1e-5).

Because each MM step reduces the surrogate, and the surrogate upper-bounds the
true loss, the algorithm guarantees that the negative log-likelihood is
non-increasing across iterations.

### New and updated functions

#### `LogBip()` -- updated

The main fitting function now accepts `method = "PDLB"` (Projection-based
logistic biplot with block coordinate Descent) in addition to the existing
`"MM"`, `"CG"`, and `"BFGS"` methods.

When `method = "PDLB"`:

* The binary matrix `x` may contain `NA` values.
* The returned object includes an `impute_x` component: the completed binary
  matrix with missing entries replaced by the model's fitted values.
* If the matrix contains missing values and a non-PDLB method is requested,
  a warning is issued and the method is automatically switched to `"PDLB"`.

```r
# Complete data -- coordinate descent MM algorithm (fast, no missing values)
res_MM <- LogBip(x = Methylation, method = "MM", maxit = 1000)

# Matrix with missing data -- projection-based block coordinate descent
set.seed(12345)
n <- nrow(Methylation); p <- ncol(Methylation)
miss         <- matrix(rbinom(n * p, 1, 0.2), n, p)
miss         <- ifelse(miss == 1, NA, miss)
x_miss       <- Methylation + miss
res_PDLB     <- LogBip(x = x_miss, method = "PDLB", maxit = 1000)
imputed_data <- res_PDLB$impute_x   # completed matrix
```

#### `proj_LogBip()` -- new

Low-level function that implements the projection-based block coordinate descent
algorithm directly. It is called internally by `LogBip(method = "PDLB")` but
is also exported for advanced users who need direct control over the algorithm.

```r
out <- proj_LogBip(x = x_miss, k = 2, max_iters = 1000, epsilon = 1e-5)
# out$mu      -- estimated intercept vector (length p)
# out$A       -- row-marker matrix (n x k)
# out$B       -- column-marker matrix (p x k)
# out$x_est   -- imputed binary matrix
# out$iter    -- number of iterations
# out$loss_funct -- loss function values per iteration
```

#### `cv_LogBip()` -- updated

Cross-validation now supports `method = "PDLB"`, allowing selection of the
optimal number of dimensions k for datasets with missing values.

```r
cv_result <- cv_LogBip(data = x_miss, k = 0:5, method = "PDLB", maxit = 1000)
```

### References

Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2021). Logistic biplot by
conjugate gradient algorithms and iterated SVD. *Mathematics*, *9*(16), 2015.
https://doi.org/10.3390/math9162015

Landgraf, A. J., & Lee, Y. (2020). Dimensionality reduction for binary data
through the projection of natural parameters. *Journal of Multivariate
Analysis*, *180*, 104668.
https://doi.org/10.1016/j.jmva.2020.104668

Pearson, K. (1901). On lines and planes of closest fit to systems of points in
space. *The London, Edinburgh, and Dublin Philosophical Magazine and Journal of
Science*, *2*(11), 559-572.

Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
M. Greenacre & J. Blasius (Eds.), *Multiple Correspondence Analysis and Related
Methods* (pp. 503-521). Chapman & Hall.

---

# BiplotML 1.0.0

## Initial CRAN release

First release of **BiplotML**, providing methods for fitting logistic biplot
models to multivariate binary data.

### Functions

* `LogBip()` -- fit a logistic biplot using conjugate gradient (`"CG"`) or
  BFGS (`"BFGS"`) optimization methods.
* `sdv_MM()` -- fit a logistic biplot via the coordinate descent
  Majorization-Minimization algorithm (`method = "MM"` in `LogBip()`).
* `bootBLB()` -- bootstrap logistic biplot with optional confidence ellipses
  for row markers.
* `plotBLB()` -- ggplot2-based biplot visualization with directed variable
  vectors and optional confidence ellipses.
* `pred_LB()` -- predict binary responses and compute per-variable optimal
  classification thresholds minimising the Balanced Error Rate.
* `fitted_LB()` -- extract fitted values on the logit or probability scale.
* `cv_LogBip()` -- k-fold cross-validation for selecting the number of
  dimensions.
* `performanceBLB()` -- compare convergence speed and accuracy across multiple
  optimization algorithms.
* `gradientDesc()` -- fit a logistic biplot via simple gradient descent
  (pedagogical / benchmarking use).
* `simBin()` -- simulate a binary data matrix from a latent variable model
  for benchmarking and cross-validation studies.

### Data

* `Methylation` -- a binary matrix of DNA methylation data (50 individuals,
  13 CpG sites) used in examples throughout the package.

### Reference

Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2021). Logistic biplot by
conjugate gradient algorithms and iterated SVD. *Mathematics*, *9*(16), 2015.
https://doi.org/10.3390/math9162015

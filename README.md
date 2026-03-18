# BiplotML

**BiplotML** implements methods for fitting logistic biplot models to multivariate binary data. A logistic biplot represents individuals as points and binary variables as directed vectors in a low-dimensional subspace. The orthogonal projection of each individual's point onto a variable's vector approximates the expected probability that the corresponding characteristic is present, providing an intuitive simultaneous visualization of observations and variables.

The package provides several fitting algorithms:

* **MM** — Coordinate descent Majorization-Minimization algorithm (fast, recommended default).
* **PDLB** — Block coordinate descent algorithm based on data projection; handles matrices with missing values and allows new individuals to be projected as supplementary rows without refitting the model.
* **CG** — Conjugate gradient algorithms (Fletcher–Reeves, Polak–Ribière, Beale–Sorenson, Dai–Yuan).
* **BFGS** — Broyden–Fletcher–Goldfarb–Shanno quasi-Newton method.

A k-fold cross-validation function (`cv\_LogBip`) is included to help select the number of dimensions.

## Installation

Install the released version from CRAN:

```r
install.packages("BiplotML")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install\_github("jgbabativam/BiplotML")
```

## Usage

```r
library(BiplotML)

data("Methylation")

# Fit a logistic biplot using the MM algorithm (default)
res\_MM <- LogBip(x = Methylation, method = "MM", maxit = 1000)

# Fit using the PDLB algorithm (supports missing data and supplementary rows)
set.seed(12345)
n <- nrow(Methylation); p <- ncol(Methylation)
miss <- matrix(rbinom(n \* p, 1, 0.2), n, p)
miss <- ifelse(miss == 1, NA, miss)
x\_miss <- Methylation + miss
res\_PDLB <- LogBip(x = x\_miss, method = "PDLB", maxit = 1000)

# Select the number of dimensions via cross-validation
cv\_result <- cv\_LogBip(data = Methylation, k = 0:5, method = "MM")

# Bootstrap confidence ellipses
set.seed(02052020)
res\_boot <- bootBLB(x = Methylation, ellipses = TRUE)
```

## Main functions

|Function|Description|
|-|-|
|`LogBip()`|Fit a logistic biplot using a chosen algorithm|
|`sdv\_MM()`|Coordinate descent MM algorithm (called internally by `LogBip`)|
|`proj\_LogBip()`|Block coordinate descent with data projection and missing-data support|
|`cv\_LogBip()`|Cross-validation to select the number of dimensions|
|`bootBLB()`|Bootstrap logistic biplot with confidence ellipses|
|`plotBLB()`|Plot a logistic biplot from a `BiplotML` object|
|`pred\_LB()`|Predict binary responses and compute optimal per-variable thresholds|
|`fitted\_LB()`|Extract fitted values on the logit or probability scale|
|`performanceBLB()`|Compare convergence and speed across multiple optimization algorithms|
|`gradientDesc()`|Fit a logistic biplot via simple gradient descent|
|`simBin()`|Simulate a binary data matrix from a latent variable model|

## Citation

If you use BiplotML in your research, please cite:

> Babativa-Márquez, J. G., \& Vicente-Villardón, J. L. (2021). Logistic biplot by
> conjugate gradient algorithms and iterated SVD. \*Mathematics\*, \*9\*(16), 2015.
> <https://doi.org/10.3390/math9162015>

## Author

Maintained by **Jose Giovany Babativa-Marquez** — jgbabativam@unal.edu.co


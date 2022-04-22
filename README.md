# *BiplotML*

Logistic Biplot is a method that allows representing multivariate binary data on a subspace of low dimension, where each individual is represented by a point and each variable as vectors directed through the origin. The orthogonal projection of individuals onto these vectors predicts the expected probability that the characteristic occurs. The package contains new techniques to estimate the model parameters and builds in each case the Logistic Biplot Model. Since version 1.1.1, an algorithm that considers missing data can be used in the LogBip function, while the cv_LogBip function implements a cross-validation algorithm that allows identifying the hyperparameter k, which represents the number of dimensions in the model. References can be found in the help for each procedure.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("jgbabativam/BiplotML")
```

## Author

This package is maintained by Giovany Babativa. Email: gbabativam@gmail.com

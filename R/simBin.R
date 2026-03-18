#' @importFrom stats qlogis
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @export
#' @title Simulate a Multivariate Binary Matrix
#' @description
#' Simulates a binary data matrix from a logistic biplot latent variable model
#' with known parameters, useful for benchmarking and cross-validation studies.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param n Number of rows (individuals).
#' @param p Number of columns (variables).
#' @param k Number of underlying latent dimensions.
#' @param D Sparsity control: the marginal probability of a 1 in the population.
#'   A value close to 0 or 1 yields a sparse or dense matrix, respectively.
#' @param C Variance scaling factor for the row scores. Default is \code{C = 1}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{X}}{Simulated binary matrix (\eqn{n \times p}).}
#'     \item{\code{P}}{Matrix of true Bernoulli probabilities (\eqn{n \times p}).}
#'     \item{\code{Theta}}{Matrix of true log-odds (natural parameters).}
#'     \item{\code{A}}{True row-marker matrix (\eqn{n \times k}).}
#'     \item{\code{B}}{True column-marker matrix (\eqn{p \times k}), orthonormal.}
#'     \item{\code{mu}}{True intercept vector of length \eqn{p}.}
#'     \item{\code{D}}{Observed proportion of ones in \strong{X}.}
#'     \item{\code{n}}{Number of rows.}
#'     \item{\code{p}}{Number of columns.}
#'   }
#'
#' @seealso \code{\link{cv_LogBip}}
#'
#' @examples
#' x <- simBin(n = 100, p = 50, k = 3, D = 0.5)

simBin <- function(n, p, k, D, C = 1) {

  for (pkg in c("pracma", "mvtnorm")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for this function. Please install it.",
           call. = FALSE)
    }
  }

  Bp <- matrix(rnorm(p * k), p)
  B  <- pracma::gramSchmidt(Bp)$Q

  S      <- diag(k)
  mu     <- rep(qlogis(D), p)
  centro <- rep(0, k)
  A      <- mvtnorm::rmvnorm(n, mean = centro, sigma = S)

  logOdds <- rep(1, n) %*% t(mu) + C * (A %*% t(B))
  P       <- plogis(logOdds)

  X <- t(sapply(seq_len(n), function(i)
    sapply(seq_len(p), function(j) rbinom(1, 1, P[i, j]))))

  desb <- sum(X) / length(X)

  list(X = X, P = P, Theta = logOdds, A = A, B = B,
       mu = mu, D = desb, n = n, p = p)
}

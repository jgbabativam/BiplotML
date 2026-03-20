#' @export
#' @title Fit a Binary Logistic Biplot with Missing Data via Block Coordinate Descent
#' @description
#' Estimates the intercept vector \eqn{\mu}, the row-marker matrix \strong{A},
#' and the column-marker matrix \strong{B} using a data-projection model with a
#' block coordinate descent algorithm. Missing values in the binary matrix are
#' imputed iteratively during model fitting. This function also allows new
#' individuals to be projected as supplementary rows without refitting the model,
#' since the row markers are derived directly from the estimated column markers.
#' This is the low-level function called by \code{\link{LogBip}} when
#' \code{method = "PDLB"}.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param x A binary matrix, possibly containing \code{NA} values.
#' @param k Number of dimensions. Default is \code{k = 2}.
#' @param max_iters Maximum number of iterations. Default is \code{1000}.
#' @param random_start Logical; if \code{TRUE}, parameters are initialised
#'   randomly. Default is \code{FALSE} (SVD initialisation).
#' @param epsilon Convergence tolerance for the relative decrease in the loss
#'   function. Default is \code{1e-5}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{mu}}{Estimated intercept vector of length \eqn{p}.}
#'     \item{\code{A}}{Estimated row-marker matrix (\eqn{n \times k}).}
#'     \item{\code{B}}{Estimated column-marker matrix (\eqn{p \times k}).}
#'     \item{\code{x_est}}{Imputed binary matrix (missing entries replaced by
#'       fitted values).}
#'     \item{\code{iter}}{Number of iterations performed.}
#'     \item{\code{loss_funct}}{Vector of normalised loss-function values at each
#'       iteration.}
#'   }
#'
#' @references
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2026). Logistic biplot
#' with missing data. \emph{In process}.
#'
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2021). Logistic biplot
#' by conjugate gradient algorithms and iterated SVD. \emph{Mathematics},
#' \emph{9}(16), 2015. \doi{10.3390/math9162015}
#'
#' Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
#' M. Greenacre & J. Blasius (Eds.), \emph{Multiple Correspondence Analysis and
#' Related Methods} (pp. 503--521). Chapman & Hall.
#'
#' @seealso \code{\link{LogBip}}, \code{\link{cv_LogBip}}
#'
#' @examples
#' \donttest{
#' data("Methylation")
#' set.seed(12345)
#' n <- nrow(Methylation); p <- ncol(Methylation)
#' miss <- matrix(rbinom(n * p, 1, 0.2), n, p)
#' miss <- ifelse(miss == 1, NA, miss)
#' x_miss <- Methylation + miss
#' out <- proj_LogBip(x = x_miss, k = 2, max_iters = 1000)
#' }

proj_LogBip <- function(x, k = 2, max_iters = 1000,
                        random_start = FALSE, epsilon = 1e-05) {

  x <- as.matrix(x)

  verify <- apply(x, 2, sd, na.rm = TRUE)
  if (any(verify == 0)) {
    stop("Some variables have zero variance; the procedure cannot be applied.")
  }

  n <- nrow(x); p <- ncol(x)

  partial_decomp <- FALSE
  if (n > 100 | p > 20) {
    partial_decomp <- TRUE
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      message("Package 'RSpectra' is recommended to speed up the algorithm.")
      partial_decomp <- FALSE
    }
  }

  if (k > p) {
    warning("k must be less than or equal to p. Setting k = p.")
    k <- p
  }

  W  <- ifelse(!is.na(x), 1, 0)
  xi <- ifelse(is.na(x), 1, x)
  x0 <- W * xi + 0.5 * (1 - W)

  if (random_start) {
    b0 <- rep(0, p)
    V  <- matrix(rnorm(p * k), p, k)
  } else {
    b0 <- colMeans(x0, na.rm = TRUE)
    if (partial_decomp) {
      udv <- RSpectra::svds(scale(x0, center = b0, scale = FALSE), k = k)
    } else {
      udv <- svd(scale(x0, center = b0, scale = FALSE))
    }
    V <- matrix(udv$v[, 1:k], p, k)
  }

  A     <- (x0 - rep(1, n) %*% t(b0)) %*% V
  B     <- V
  theta <- rep(1, n) %*% t(b0) + A %*% t(B)
  P     <- plogis(theta)

  thresh <- thresholds(x0, P, ncuts = 100)$thres

  err    <- numeric(max_iters + 1)
  xl     <- x0
  err[1] <- log_like(x = xl, w = W, theta = theta) / (n * p)

  l <- 0
  for (l in seq_len(max_iters)) {
    Zl  <- theta + 4 * (xl - P)
    Ml  <- W * Zl + (1 - W) * theta
    Xlc <- xl - rep(1, n) %*% t(b0)
    Mlc <- Ml - rep(1, n) %*% t(b0)
    Yl  <- -(t(Xlc) %*% Xlc - t(Xlc) %*% Mlc - t(Mlc) %*% Xlc)

    sdv <- RSpectra::eigs_sym(Yl, k, which = "LM")
    if (any(sdv$values[1:k] < 0)) {
      sdv <- eigen(Yl, symmetric = TRUE)
    }

    V  <- matrix(sdv$vectors[, 1:k], p, k)
    b0 <- as.vector(1 / n * t(Ml - xl %*% tcrossprod(V)) %*% rep(1, n))

    theta <- rep(1, n) %*% t(b0) +
      (xl - rep(1, n) %*% t(b0)) %*% tcrossprod(V)
    P <- plogis(theta)

    alphas <- rep(1, n) %*% t(as.matrix(thresh$threshold))
    xhat   <- ifelse(P < alphas, 0, 1)
    xl     <- W * xi + (1 - W) * xhat

    err[l + 1] <- log_like(x = xl, w = W, theta = theta) /
      (nrow(x) * ncol(x))

    if (l == max_iters) {
      warning("The algorithm reached ", max_iters,
              " iterations without converging.")
    }

    if (l > 2) {
      if ((err[l] - err[l + 1]) / err[l] < epsilon) break
    }
  }

  A <- (xl - rep(1, n) %*% t(b0)) %*% V
  B <- V

  out <- list(mu = b0, A = A, B = B, x_est = xl,
              iter = l, loss_funct = err[seq_len(l)])
  return(out)
}

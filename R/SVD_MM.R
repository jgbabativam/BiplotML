#' @export
#' @title Fit a Binary Logistic Biplot via Coordinate Descent MM Algorithm
#' @description
#' Estimates the intercept vector \eqn{\mu}, the row-marker matrix \strong{A},
#' and the column-marker matrix \strong{B} using an iterative coordinate descent
#' Majorization-Minimization (MM) algorithm. This is the low-level function
#' called by \code{\link{LogBip}} when \code{method = "MM"}.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param x A binary matrix with no missing values.
#' @param k Number of dimensions. Default is \code{k = 5}.
#' @param iterations Maximum number of iterations. Default is \code{1000}.
#' @param truncated Logical; if \code{TRUE} (default for large matrices), the
#'   truncated SVD from \pkg{RSpectra} is used to speed up computation.
#' @param random Logical; if \code{TRUE}, parameters are initialised randomly.
#'   Default is \code{FALSE} (SVD initialisation).
#' @param epsilon Convergence tolerance. The algorithm stops when the relative
#'   decrease in the loss function is below this value. Default is \code{1e-4}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{mu}}{Estimated intercept vector of length \eqn{p}.}
#'     \item{\code{A}}{Estimated row-marker matrix (\eqn{n \times k}).}
#'     \item{\code{B}}{Estimated column-marker matrix (\eqn{p \times k}).}
#'     \item{\code{iterations}}{Number of iterations performed.}
#'     \item{\code{loss_func}}{Vector of normalised loss-function values at each
#'       iteration.}
#'   }
#'
#' @references
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
#' out <- sdv_MM(x = Methylation)
#' }
sdv_MM <- function(x, k = 5, iterations = 1000, truncated = TRUE,
                   random = FALSE, epsilon = 1e-4) {

  n <- nrow(x); p <- ncol(x)

  if (k > p) {
    message("k must be less than or equal to the number of columns. Setting k = p.")
    k <- p
  }

  if (!random) {
    mu <- colMeans(8 * x)
    if (truncated) {
      vp <- RSpectra::svds(scale(8 * x, center = TRUE, scale = FALSE),
                           k = k, nu = k, nv = k)
    } else {
      vp <- svd(scale(8 * x, center = TRUE, scale = FALSE))
    }
    if (k == 1) {
      A <- matrix(vp$u[, 1:k], n, k) * vp$d[1]
      B <- matrix(vp$v[, 1:k], p, k)
    } else {
      A <- vp$u[, 1:k] %*% diag(vp$d[1:k])
      B <- vp$v[, 1:k]
    }
  } else {
    mu <- runif(p)
    A  <- matrix(runif(n * k), n, k)
    B  <- matrix(runif(p * k), p, k)
  }

  theta <- rep(1, n) %*% t(mu) + (A %*% t(B))
  loss_func <- numeric(iterations + 1)
  pi_hat <- plogis(theta)
  loss_func[1] <- -sum(x * log(pi_hat) + (1 - x) * log(1 - pi_hat),
                       na.rm = TRUE) / (n * p)

  j <- 0
  for (i in seq_len(iterations)) {
    old_mu <- mu; old_A <- A; old_B <- B

    Z  <- theta + 4 * (x - pi_hat)
    mu <- colMeans(Z)

    if (truncated) {
      vp <- RSpectra::svds(scale(Z, center = TRUE, scale = FALSE),
                           k = k, nu = k, nv = k)
    } else {
      vp <- svd(scale(Z, center = TRUE, scale = FALSE))
    }

    if (k == 1) {
      A <- matrix(vp$u[, 1:k], n, k) * vp$d[1]
      B <- matrix(vp$v[, 1:k], p, k)
    } else {
      A <- vp$u[, 1:k] %*% diag(vp$d[1:k])
      B <- vp$v[, 1:k]
    }

    theta   <- rep(1, n) %*% t(mu) + (A %*% t(B))
    pi_hat  <- plogis(theta)
    loss_func[i + 1] <- -sum(x * log(pi_hat) + (1 - x) * log(1 - pi_hat),
                             na.rm = TRUE) / (n * p)

    if (loss_func[i] < loss_func[i + 1]) {
      mu <- old_mu; A <- old_A; B <- old_B
      if (truncated) {
        warning("The loss function increased at iteration ", i,
                ". Consider rerunning with truncated = FALSE.")
      }
    }

    if (i > 10) {
      if ((loss_func[i] - loss_func[i + 1]) / loss_func[i] < epsilon)
        break
    }

    if (i == iterations) {
      warning("The algorithm reached the maximum of ", iterations,
              " iterations without converging.")
    }

    j <- i
  }

  out <- list(mu = mu, A = A, B = B, iterations = j,
              loss_func = loss_func[seq_len(j)])
  return(out)
}

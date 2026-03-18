#' @export
#' @title Fit a Binary Logistic Biplot via Gradient Descent
#' @description
#' Estimates the row-marker matrix \strong{A} and the column-marker matrix
#' \strong{B} of a binary logistic biplot using a simple (batch) gradient
#' descent algorithm. This function is mainly provided for pedagogical purposes
#' and benchmarking; the MM and CG methods in \code{\link{LogBip}} are
#' generally faster and more reliable.
#'
#' @details
#' The model is
#' \deqn{\mathrm{logit}(\pi_{ij}) =
#'   \log\!\left(\frac{\pi_{ij}}{1-\pi_{ij}}\right) =
#'   \mu_j + \sum_{s=1}^k b_{js}\,a_{is} = \mu_j + \mathbf{a}_i^\top \mathbf{b}_j.}
#' The gradient with respect to the full parameter vector is
#' \deqn{\nabla\ell =
#'   \left(\frac{\partial\ell}{\partial\boldsymbol{\mu}},\,
#'         \frac{\partial\ell}{\partial\mathbf{A}},\,
#'         \frac{\partial\ell}{\partial\mathbf{B}}\right) =
#'   \left((\boldsymbol{\Pi}-\mathbf{X})^\top,\;
#'         (\boldsymbol{\Pi}-\mathbf{X})\mathbf{B},\;
#'         (\boldsymbol{\Pi}-\mathbf{X})^\top\mathbf{A}\right).}
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param x A binary matrix.
#' @param k Number of dimensions. Default is \code{k = 2}.
#' @param rate Learning rate \eqn{\alpha} for the gradient descent update.
#'   Default is \code{0.001}.
#' @param converg Convergence tolerance: the algorithm stops when the relative
#'   change in the loss function is below this value. Default is \code{0.001}.
#' @param max_iter Maximum number of iterations.
#' @param plot Logical; if \code{TRUE}, the logistic biplot is plotted after
#'   fitting. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{BiplotML} (a named list) containing:
#'   \describe{
#'     \item{\code{Ahat}}{Estimated row-marker matrix.}
#'     \item{\code{Bhat}}{Estimated column-marker matrix (including intercepts).}
#'     \item{\code{method}}{Character string \code{"Gradient Descent"}.}
#'   }
#'
#' @references
#' Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
#' M. Greenacre & J. Blasius (Eds.), \emph{Multiple Correspondence Analysis and
#' Related Methods} (pp. 503--521). Chapman & Hall.
#'
#' @seealso \code{\link{plotBLB}}, \code{\link{performanceBLB}}
#'
#' @examples
#' data("Methylation")
#' set.seed(02052020)
#' outGD <- gradientDesc(x = Methylation, k = 2, max_iter = 10000, plot = TRUE)

gradientDesc <- function(x, k = 2, rate = 0.001, converg = 0.001,
                         max_iter, plot = FALSE, ...) {

  x   <- as.matrix(x)
  n   <- nrow(x); p <- ncol(x)
  aik <- n * k;   bjk <- p * (k + 1)
  dTheta <- aik + bjk

  par <- runif(dTheta)
  A   <- matrix(par[1:aik], n, k)
  B   <- matrix(par[(aik + 1):dTheta], p, k + 1)

  lin  <- cbind(rep(1, n), A) %*% t(B)
  hX   <- plogis(lin)
  J    <- -sum(x * log(hX) + (1 - x) * log(1 - hX), na.rm = TRUE)
  MSE  <- J / n

  converged  <- FALSE
  iterations <- 0

  while (!converged) {
    A_new <- A - rate * (hX - x) %*% B[, -1]
    B_new <- B - rate * t(hX - x) %*% cbind(rep(1, n), A)

    lin_new <- cbind(rep(1, n), A_new) %*% t(B_new)
    hX_hat  <- plogis(lin_new)
    J_new   <- -sum(x * log(hX_hat) + (1 - x) * log(1 - hX_hat),
                    na.rm = TRUE)
    MSE_new <- J_new / n

    iterations <- iterations + 1

    if (abs(MSE_new - MSE) <= converg) {
      converged <- TRUE
    }

    if (iterations >= max_iter) {
      warning("The algorithm reached ", max_iter,
              " iterations without converging.")
      converged <- TRUE
    }

    A   <- A_new; B <- B_new
    hX  <- hX_hat
    MSE <- MSE_new
  }

  Ahat <- data.frame(A)
  Bhat <- data.frame(B)
  colnames(Ahat) <- paste0("Dim", seq_len(k))
  colnames(Bhat) <- paste0("bb",  seq(0, k, 1))
  rownames(Ahat) <- rownames(x)
  rownames(Bhat) <- colnames(x)

  out <- list(Ahat = Ahat, Bhat = Bhat, method = "Gradient Descent")

  if (plot & ncol(Ahat) > 1) {
    print(plotBLB(x = out, ellipses = FALSE,
                  titles = "Logistic Biplot - Gradient Descent"))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

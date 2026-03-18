#' @importFrom optimx optimr
#' @importFrom stats runif
#' @importFrom stats sd
#' @export
#' @title Fit a Binary Logistic Biplot
#' @description
#' Estimates the intercept vector \eqn{\mu}, the row-marker matrix \strong{A}, and the
#' column-marker matrix \strong{B} of a logistic biplot model using the optimization
#' algorithm selected by the user.
#'
#' @details
#' The following fitting methods are available:
#'
#' \strong{Conjugate gradient (CG):} Set \code{method = "CG"} and choose the update
#' formula via \code{type}:
#' \itemize{
#'   \item \code{type = 1} --- Fletcher--Reeves
#'   \item \code{type = 2} --- Polak--Ribiere
#'   \item \code{type = 3} --- Hestenes--Stiefel
#'   \item \code{type = 4} --- Dai--Yuan
#' }
#'
#' \strong{Coordinate descent MM:} Set \code{method = "MM"} to use the iterative
#' coordinate descent Majorization-Minimization algorithm.
#'
#' \strong{Projection-based algorithm (PDLB):} Set \code{method = "PDLB"} when the
#' binary matrix contains missing values, or when the row coordinates of new
#' (supplementary) individuals need to be estimated without refitting the model.
#' See Babativa-Marquez & Vicente-Villardon (2022) for details.
#'
#' \strong{BFGS:} Set \code{method = "BFGS"} to use the Broyden--Fletcher--
#' Goldfarb--Shanno quasi-Newton method.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param x A binary matrix (or a matrix with \code{NA} values when
#'   \code{method = "PDLB"}).
#' @param k Number of dimensions. Default is \code{k = 2}.
#' @param method Fitting algorithm. One of \code{"MM"} (default), \code{"CG"},
#'   \code{"PDLB"}, or \code{"BFGS"}.
#' @param type Update formula for the conjugate gradient method: \code{1} =
#'   Fletcher--Reeves, \code{2} = Polak--Ribiere, \code{3} = Beale--Sorenson.
#'   Ignored for other methods.
#' @param plot Logical; if \code{TRUE} (default), the logistic biplot is plotted
#'   after fitting.
#' @param maxit Maximum number of iterations. Defaults to \code{100} for gradient
#'   methods and \code{500} for derivative-free methods.
#' @param endsegm End point of the variable segment on the probability scale.
#'   The segment starts at 0.5 and ends at this value. Default is \code{0.90}.
#' @param label.ind Logical; if \code{TRUE}, row points are labelled. Default is
#'   \code{FALSE}.
#' @param col.ind Color for the row markers. Passed to \code{\link{plotBLB}}.
#' @param draw Which graph to draw: \code{"biplot"} (default) for both rows and
#'   columns, \code{"ind"} for individuals only, or \code{"var"} for variables
#'   only.
#' @param random_start Logical; if \code{TRUE}, parameters are initialised
#'   randomly. If \code{FALSE} (default), an SVD-based initialisation is used.
#' @param L Ridge penalization parameter. Default is \code{L = 0} (no penalty).
#' @param cv_LogBip Logical; indicates whether the function is being called
#'   internally by \code{\link{cv_LogBip}}. Users should leave this as
#'   \code{FALSE} (default).
#'
#' @return An object of class \code{BiplotML} (a named list) containing:
#'   \describe{
#'     \item{\code{Ahat}}{Data frame of row-marker coordinates.}
#'     \item{\code{Bhat}}{Data frame of column-marker coordinates, including the
#'       intercept column \code{bb0}.}
#'     \item{\code{method}}{Character string identifying the fitting method used.}
#'     \item{\code{loss_function}}{Vector of loss-function values at each
#'       iteration (MM and PDLB methods only).}
#'     \item{\code{iterations}}{Number of iterations performed (MM and PDLB
#'       methods only).}
#'     \item{\code{impute_x}}{Imputed binary matrix (PDLB method only).}
#'   }
#'
#' @references
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2022). Logistic biplot
#' with missing data. \emph{In process}.
#'
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2021). Logistic biplot
#' by conjugate gradient algorithms and iterated SVD. \emph{Mathematics},
#' \emph{9}(16), 2015. \doi{10.3390/math9162015}
#'
#' Nash, J. C. (2011). Unifying optimization algorithms to aid software system
#' users: optimx for R. \emph{Journal of Statistical Software}, \emph{43}(9),
#' 1--14.
#'
#' Nash, J. C. (2014). On best practice optimization methods in R.
#' \emph{Journal of Statistical Software}, \emph{60}(2), 1--14.
#'
#' Nocedal, J., & Wright, S. (2006). \emph{Numerical Optimization} (2nd ed.).
#' Springer.
#'
#' Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
#' M. Greenacre & J. Blasius (Eds.), \emph{Multiple Correspondence Analysis and
#' Related Methods} (pp. 503--521). Chapman & Hall.
#'
#' @seealso \code{\link{plotBLB}}, \code{\link{pred_LB}}, \code{\link{fitted_LB}}
#'
#' @examples
#' \donttest{
#' data("Methylation")
#'
#' # Fit using the coordinate descent MM algorithm
#' res_MM <- LogBip(x = Methylation, method = "MM", maxit = 1000)
#'
#' # Fit using the PDLB algorithm with simulated missing data
#' set.seed(12345)
#' n <- nrow(Methylation); p <- ncol(Methylation)
#' miss <- matrix(rbinom(n * p, 1, 0.2), n, p)
#' miss <- ifelse(miss == 1, NA, miss)
#' x_miss <- Methylation + miss
#' res_PDLB <- LogBip(x = x_miss, method = "PDLB", maxit = 1000)
#' }

LogBip <- function(x, k = 2, method = "MM", type = NULL, plot = TRUE,
                   maxit = NULL, endsegm = 0.90, label.ind = FALSE,
                   col.ind = NULL,
                   draw = c("biplot", "ind", "var"),
                   random_start = FALSE, L = 0, cv_LogBip = FALSE) {

  x <- as.matrix(x)
  verify <- apply(x, 2, sd, na.rm = TRUE)

  if (any(verify == 0)) {
    stop("Some variables have zero variance; the procedure cannot be applied.")
  }

  if (cv_LogBip & any(is.na(x))) {
    x <- sweep(x, MARGIN = 2,
               STATS = colMeans(x, na.rm = TRUE),
               FUN = function(z, s) ifelse(is.na(z), s, z))
  } else if (any(is.na(x)) & method != "PDLB") {
    warning("The binary matrix contains missing values; method '", method,
            "' has been switched to 'PDLB'.")
    method <- "PDLB"
  }

  n <- nrow(x); p <- ncol(x); aik <- n * k; bjk <- p * (k + 1)
  dTheta <- aik + bjk; s <- k + 1

  truncated <- FALSE
  if (n > 100 | p > 20) {
    truncated <- TRUE
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      message("Package 'RSpectra' is recommended to speed up the algorithm.")
      truncated <- FALSE
    }
  }

  if (method != "PDLB") {
    if (random_start) {
      params <- runif(dTheta)
    } else {
      b0 <- colMeans(4 * x, na.rm = TRUE)
      if (truncated) {
        udv <- RSpectra::svds(scale(4 * x, center = b0, scale = FALSE),
                              k = k)
      } else {
        udv <- svd(scale(4 * x, center = b0, scale = FALSE))
      }
      A <- matrix(udv$u[, 1:k], n, k) %*% diag(udv$d[1:k], nrow = k, ncol = k)
      B <- matrix(udv$v[, 1:k], p, k)
      params <- c(A, b0, B)
    }
  }

  if (method == "CG") {
    res <- .run_optimr(par = params, fn = J.BipLog.BIN,
                       gr = Grad.BipLog.BIN,
                       xt = x, k = k, lambda = L,
                       method = method,
                       ctrl = list(type = type))
  } else if (method == "MM") {
    if (!is.null(maxit)) {
      res <- sdv_MM(x = x, k = k, iterations = maxit,
                    random = random_start, truncated = truncated)
    } else {
      res <- sdv_MM(x = x, k = k, random = random_start,
                    truncated = truncated)
    }
  } else if (method == "PDLB") {
    if (!is.null(maxit)) {
      res <- proj_LogBip(x = x, k = k, max_iters = maxit,
                         random_start = random_start)
    } else {
      res <- proj_LogBip(x = x, k = k, random_start = random_start)
    }
  } else {
    res <- .run_optimr(par = params, fn = J.BipLog.BIN,
                       gr = Grad.BipLog.BIN,
                       xt = x, k = k, lambda = L,
                       method = method)
  }

  if (method %in% c("MM", "PDLB")) {
    Ahat <- data.frame(res$A)
  } else {
    par <- res$par
    Ahat <- data.frame(matrix(res$par[1:aik], n, k))
  }
  colnames(Ahat) <- paste0("Dim", seq_len(k))
  rownames(Ahat) <- rownames(x)

  if (method %in% c("MM", "PDLB")) {
    Bhat <- data.frame(res$mu, res$B)
  } else {
    Bhat <- data.frame(matrix(par[(aik + 1):dTheta], p, k + 1))
  }
  rownames(Bhat) <- colnames(x)
  colnames(Bhat) <- paste0("bb", seq(0, k, 1))

  if (method == "MM") {
    method <- "coordinate descent MM"
    out <- list(Ahat = Ahat, Bhat = Bhat, method = method,
                loss_function = res$loss_func, iterations = res$iterations)
  } else if (method == "PDLB") {
    method <- "data projection with block coordinate descent"
    out <- list(Ahat = Ahat, Bhat = Bhat, impute_x = res$x_est,
                method = method,
                loss_function = res$loss_funct, iterations = res$iter)
  } else {
    if (method == "CG" & type == 1) method <- "CG: Fletcher--Reeves"
    if (method == "CG" & type == 2) method <- "CG: Polak--Ribiere"
    if (method == "CG" & type == 3) method <- "CG: Beale--Sorenson"
    if (method == "Rcgmin")         method <- "CG: Dai--Yuan"
    out <- list(Ahat = Ahat, Bhat = Bhat, method = method)
  }

  if (plot & ncol(Ahat) > 1) {
    print(plotBLB(x = out, ellipses = FALSE, endsegm = endsegm,
                  titles = "Logistic Biplot", label.ind = label.ind,
                  col.ind = col.ind, draw = draw))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

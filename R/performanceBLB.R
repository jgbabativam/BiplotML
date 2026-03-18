#' @importFrom optimx optimx
#' @export
#' @title Compare Optimization Algorithms for Binary Logistic Biplot Estimation
#' @description
#' Fits the binary logistic biplot model using multiple optimization algorithms
#' and returns a summary of their computation time, convergence status, and
#' number of function evaluations, facilitating algorithm selection.
#'
#' @details
#' The following algorithm groups are available via the \code{method} argument:
#' \itemize{
#'   \item \code{1} --- Derivative-free methods: Nelder-Mead, UOBYQA, NEWUOA.
#'   \item \code{2} --- Gradient methods (default): CG, Rcgmin.
#'   \item \code{3} --- Quasi-Newton methods: BFGS, L-BFGS-B, nlm, nlminb.
#'   \item \code{4} --- All of the above.
#' }
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param xi A binary matrix.
#' @param k Number of dimensions. Default is \code{k = 2}.
#' @param L Ridge penalization parameter. Default is \code{L = 0}.
#' @param method Algorithm group to compare: \code{1} (derivative-free),
#'   \code{2} (gradient, default), \code{3} (quasi-Newton), or \code{4} (all).
#' @param maxit Maximum number of iterations per algorithm.
#'
#' @return A data frame with one row per algorithm and columns:
#'   \describe{
#'     \item{\code{method}}{Algorithm name.}
#'     \item{\code{evaluat}}{Final value of the objective function.}
#'     \item{\code{convergence}}{Convergence status.}
#'     \item{\code{fevals}}{Number of function evaluations.}
#'     \item{\code{time}}{Elapsed computation time.}
#'   }
#'
#' @references
#' Nash, J. C. (2011). Unifying optimization algorithms to aid software system
#' users: optimx for R. \emph{Journal of Statistical Software}, \emph{43}(9),
#' 1--14.
#'
#' Nash, J. C. (2014). On best practice optimization methods in R.
#' \emph{Journal of Statistical Software}, \emph{60}(2), 1--14.
#'
#' Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
#' M. Greenacre & J. Blasius (Eds.), \emph{Multiple Correspondence Analysis and
#' Related Methods} (pp. 503--521). Chapman & Hall.
#'
#' @seealso \code{\link{gradientDesc}}
#'
#' @examples
#' \donttest{
#' data("Methylation")
#' set.seed(123456)
#'
#' # Gradient methods (default)
#' performanceBLB(xi = Methylation)
#' performanceBLB(xi = Methylation, maxit = 150)
#'
#' # Derivative-free methods
#' performanceBLB(xi = Methylation, method = 1)
#' performanceBLB(xi = Methylation, method = 1, maxit = 100)
#'
#' # Quasi-Newton methods
#' performanceBLB(xi = Methylation, method = 3)
#' performanceBLB(xi = Methylation, method = 3, maxit = 100)
#'
#' # All methods
#' performanceBLB(xi = Methylation, method = 4)
#' }

performanceBLB <- function(xi, k = 2, L = 0, method = NULL, maxit = NULL) {

  n      <- nrow(xi)
  p      <- ncol(xi)
  aik    <- n * k
  bjk    <- p * (k + 1)
  dTheta <- aik + bjk

  if (is.null(method)) method <- 2

  method0 <- switch(as.character(method),
    "1" = c("Nelder-Mead", "uobyqa", "newuoa"),
    "2" = c("CG", "Rcgmin"),
    "3" = c("BFGS", "L-BFGS-B", "nlm", "nlminb"),
    "4" = c("Nelder-Mead", "CG", "BFGS", "nlm", "nlminb",
             "Rcgmin", "uobyqa", "newuoa"),
    stop("'method' must be 1, 2, 3, or 4.")
  )

  args <- list(par = runif(dTheta), fn = J.BipLog.BIN,
               gr = Grad.BipLog.BIN, xt = xi, k = k,
               lambda = L, method = method0)
  if (!is.null(maxit)) args$itnmax <- maxit

  old_wd <- setwd(tempdir())
  on.exit({
    for (f in c(file.path(tempdir(), "badhess.txt"),
                file.path(old_wd,    "badhess.txt"))) {
      if (file.exists(f)) unlink(f, force = TRUE)
    }
    setwd(old_wd)
  }, add = TRUE)
  res <- suppressWarnings(do.call(optimx, args))

  res$result <- ifelse(
    res$convcode == 0,    "convergence",
    ifelse(res$convcode == 1, "maximum iterations reached",
    ifelse(res$convcode %in% c(10, 20, 21), "inadmissible value",
    ifelse(res$convcode == 9999, "method failed", NA))))

  out <- data.frame(
    method      = rownames(res),
    evaluat     = res$value,
    convergence = res$result,
    fevals      = res$fevals,
    time        = res$xtime,
    stringsAsFactors = FALSE
  )
  return(out)
}

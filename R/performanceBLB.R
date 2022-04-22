#' @importFrom optimx optimx
#' @export
#'
#' @title
#' Performance comparison of severals estimation algorithms
#' @description
#' This function computes the estimates of A and B matrix with severals algorithms.
#' @return
#' data frame with method, time of process, convergence and number of evaluations
#' @details
#' This function compare the process time and convergence of different algorithms without gradient, with gradient or quasi-newton method for estimating the parameters in a Binary Logistic Biplot
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param xi Binary matrix.
#' @param k Dimensions number. By default \code{k = 2}.
#' @param L Penalization parameter. By default \code{L = 0}.
#' @param method use value 1 for algorithms without gradient, 2 with gradient, 3 quasi-newton methods or 4 for all methods. By default \code{method = 2}.
#' @param maxit The maximum number of iterations. Defaults to 100 for the gradient methods, and 500 without gradient.
#' @references
#' John C. Nash (2011). Unifying Optimization Algorithms to Aid Software System Users:optimx for R. Journal of Statistical Software. 43(9). 1--14.
#'
#' John C. Nash (2014). On Best Practice Optimization Methods in R. Journal of Statistical Software. 60(2). 1--14.
#'
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{gradientDesc}}
#' @examples
#' \donttest{
#' data('Methylation')
#' set.seed(123456)
#' ########### Gradient Methods
#' performanceBLB(xi = Methylation)
#' performanceBLB(xi = Methylation, maxit = 150)
#'
#' ########### Without Gradient Methods
#' performanceBLB(xi = Methylation, method = 1)
#' performanceBLB(xi = Methylation, method = 1, maxit = 100)
#'
#' ########### Quasi-Newton Methods
#' performanceBLB(xi = Methylation, method = 3)
#' performanceBLB(xi = Methylation, method = 3, maxit = 100)
#'
#' ########### All methods
#' performanceBLB(x = Methylation, method = 4)
#' }

performanceBLB <- function(xi, k = 2, L = 0, method = NULL, maxit = NULL) {
    n = nrow(xi)
    p = ncol(xi)
    aik = n * k
    bjk = p * (k + 1)
    dTheta = aik + bjk
    if (is.null(method))
        method = 2

    if (method == 1 & !is.null(maxit)) {
        method0 = c("Nelder-Mead", "uobyqa", "newuoa")
        res = optimx(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
            xt = xi, k = k, lambda = L, method = method0, itnmax = maxit)
    }
    if (method == 1 & is.null(maxit)) {
        method0 = c("Nelder-Mead", "uobyqa", "newuoa")
        res = optimx(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
            xt = xi, k = k, lambda = L, method = method0)
    }

    if (method == 2 & !is.null(maxit)) {
        method0 <- c("CG", "Rcgmin")
        res <- optimx(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
            xt = xi, k = k, lambda = L, method = method0, itnmax = maxit)

    }
    if (method == 2 & is.null(maxit)) {
        method0 = c("CG", "Rcgmin")
        res = optimx(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
            xt = xi, k = k, lambda = L, method = method0)
    }

    if (method == 3 & !is.null(maxit)) {
        method0 = c("BFGS", "L-BFGS-B", "nlm", "nlminb")
        res = optimx(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
            xt = xi, k = k, lambda = L, method = method0, itnmax = maxit)
    }
    if (method == 3 & is.null(maxit)) {
        method0 = c("BFGS", "L-BFGS-B", "nlm", "nlminb")
        res = optimx(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
            xt = xi, k = k, lambda = L, method = method0)
    }
    if (method == 4) {
        method0 = c("Nelder-Mead", "CG", "BFGS", "nlm", "nlminb", "Rcgmin", "uobyqa", "newuoa")
        res = optimx(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
            xt = xi, k = k, lambda = L, method = method0)
    }
    res$result = ifelse(res$convcode == 0, "convergence", ifelse(res$convcode ==
        1, "max iterations", ifelse(res$convcode %in% c(10, 20, 21), "inadmissible",
        ifelse(res$convcode == 9999, "method has failed", NA))))
    out = data.frame(method = rownames(res), evaluat = res$value, convergence = res$result,
        fevals = res$fevals, time = res$xtime)

    return(out)
}

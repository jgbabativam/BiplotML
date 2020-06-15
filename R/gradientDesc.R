#' @export
#' @title
#' Gradient function for Binary Logistic Biplot
#' @description
#' This function computes the parameters of A and B in Binary Logistic Biplot under algorithm of Descendent Gradient.
#' @return
#' The coefficients of A and B matrix.
#' @details
#' We note that the Binary Logistic Biplot is defined as: \deqn{logit(\pi_{ij})= log\left( \frac{\pi_{ij}}{1-\pi_{ij}} \right)=\mu_{j}+\sum_{s=1}^kb_{js}a_{is} = \mu_{j}+\mathbf{a_i^{T}b_j}}
#' Also, note that the gradient is: \deqn{\nabla \ell= \left(\frac{\partial \ell}{\partial \mu}, \frac{\partial \ell}{\partial \mathbf{A}}, \frac{\partial \ell}{\partial \mathbf{B}}\right)== \left( (\Pi - \mathbf{X})^T,  (\Pi - \mathbf{X})\mathbf{B}, (\Pi - \mathbf{X})^TA  \right)}
#'
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x Binary matrix.
#' @param k Dimensions number. By default \code{k = 2}.
#' @param rate The value of the rate of descent \eqn{\alpha} in the algorithm of descending gradient. By default \eqn{\alpha = 0.001}.
#' @param converg Tolerance limit to achieve convergence. By default \code{converg = 0.001}
#' @param max_iter Maximum iterations number.
#' @param plot Plot the Logistic Biplot.
#' @param ... other arguments
#' @references
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{plotBLB}, \link{performanceBLB}}
#' @examples
#' data('Methylation')
#' set.seed(02052020)
#' MatGD <- gradientDesc(x = Methylation, k=2, max_iter=10000)
#' outGD <- gradientDesc(x = Methylation, k=2, max_iter=10000, plot = TRUE)

gradientDesc <- function(x, k = 2, rate = 0.001, converg = 0.001, max_iter,
                         plot = FALSE, ...) {

    x = as.matrix(x)
    n = nrow(x)
    p = ncol(x)
    aik = n * k
    bjk = p * (k + 1)

    dTheta = aik + bjk
    par = runif(dTheta)
    A = matrix(par[1:aik], n, k)
    B = matrix(par[(aik + 1):dTheta], p, k + 1)

    lin = cbind(rep(1, n), A) %*% t(B)
    hX = exp(lin)/(1 + exp(lin))

    J = -sum(x * log(hX) + (1 - x) * log(1 - hX), na.rm = TRUE)

    MSE <- J/n
    converged = F
    iterations = 0
    while (converged == F) {
        A_new <- A - rate * (hX - x) %*% B[, -1]
        B_new <- B - rate * t(hX - x) %*% cbind(rep(1, n), A)

        lin_new = cbind(rep(1, n), A_new) %*% t(B_new)
        hX_hat = exp(lin_new)/(1 + exp(lin_new))

        J_new = -sum(x * log(hX_hat) + (1 - x) * log(1 - hX_hat), na.rm = TRUE)
        MSE_new <- J_new/n

        A <- A_new
        B <- B_new
        if (MSE - MSE_new <= converg) {
            converged = T
            A <- as.data.frame(A)
            B <- as.data.frame(B)
            colnames(B) = c(paste0("bb", seq(0,k,1)))
            rownames(A) = rownames(x)
            colnames(A) = c(paste0("Dim", seq(1,k,1)))
            out = list(Ahat = A, Bhat = B, method="Gradient descent")
            print(paste("The process converge with", iterations, "iterations"))

            if (plot & ncol(A)>1) {
                print(plotBLB(x=out))
            }

            class(out) <- c("BiplotML", "list")
            return(out)
        }
        iterations = iterations + 1
        if (iterations > max_iter) {
            converged = T
            A <- as.data.frame(A)
            B <- as.data.frame(B)
            colnames(B) = c(paste0("bb", seq(0,k,1)))
            rownames(A) = rownames(x)
            colnames(A) = c(paste0("Dim", seq(1,k,1)))
            out = list(Ahat = A, Bhat = B, method="Gradient descent")
            print(paste("The process not converge with", max_iter, "iterations"))

            if (plot & ncol(A)>1) {
                print(plotBLB(x=out))
            }
            class(out) <- c("BiplotML", "list")
            return(out)

        }
    }

}


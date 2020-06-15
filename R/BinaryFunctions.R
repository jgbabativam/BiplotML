#' @export
#'
#' @title
#' Necessary functions
#' @description
#' These functions are necessary for the correct performance of the package.
#' @return
#' Output functions
#' @details
#' These functions are necessary for the correct performance of the package.
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param xt Input matrix.
#' @param par Parameters.
#' @param k Dimensions number.
#' @param lambda Penalization parameter.
#' @param ... other arguments
#' @export

J.BipLog.BIN <- function(xt, par, k, lambda, ...) {
    #---- Para iniciar.
    x = as.matrix(xt)
    n = nrow(x)
    p = ncol(x)
    aik = n * k
    bjk = p * (k + 1)
    dTheta = aik + bjk

    #--Matrices de parámetros.

    A = matrix(par[1:aik], n, k)
    B = matrix(par[(aik + 1):dTheta], p, k + 1)

    #--- Parte lineal.
    lin = cbind(rep(1, n), A) %*% t(B)

    #--- Hipotesis.
    hX = exp(lin)/(1 + exp(lin))

    #---- Minimizar:
    J = -sum(x * log(hX) + (1 - x) * log(1 - hX), na.rm = TRUE) + lambda *
        sum(A^2, na.rm = TRUE)/2 + lambda * sum(B[, -1]^2, na.rm = TRUE)/2
    return(J)
}

#' @title
#' Necessary functions for gradients
#' @description
#' These functions are necessary for the correct performance of the package.
#' @return
#' Output functions
#' @details
#' These functions are necessary for the correct performance of the package.
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param xt Input matrix.
#' @param par Parameters.
#' @param k Dimensions number.
#' @param lambda Penalization parameter.
#' @param ... other arguments
#' @export

Grad.BipLog.BIN <- function(xt, par, k, lambda, ...) {
    x = as.matrix(xt)
    n = nrow(x)
    p = ncol(x)
    aik = n * k
    bjk = p * (k + 1)
    dTheta = aik + bjk

    A = matrix(par[1:aik], n, k)
    B = matrix(par[(aik + 1):dTheta], p, k + 1)

    lin = cbind(rep(1, n), A) %*% t(B)
    hX = exp(lin)/(1 + exp(lin))

    dA = (hX - x) %*% B[, -1] + lambda * A
    dB = t(hX - x) %*% cbind(rep(1, n), A) + lambda * cbind(rep(0, p), B[,
        -1])
    gradient = c(c(dA), c(dB))
    return(gradient)
}

#' @title
#' Supplementary idividuals
#' @param xs Input matrix supplematary individuals.
#' @param B Resamples number.
#' @param k Dimensions number.
#' @param par Parameters.
#' @param lambda Penalization parameter.
#' @param ... other arguments

#' @export
Indsup.BIN <- function(xs, B, par, k, lambda, ...) {
    x = as.matrix(xs)
    n = nrow(x)
    p = ncol(x)
    aik = n * k
    dTheta = aik

    A = matrix(par[1:aik], n, k)
    Bs = as.matrix(B)
    lin = cbind(rep(1, n), A) %*% t(Bs)

    hX = exp(lin)/(1 + exp(lin))

    L = -sum(x * log(hX) + (1 - x) * log(1 - hX), na.rm = TRUE) + lambda *
        sum(A^2, na.rm = TRUE)/2 + lambda * sum(Bs[, -1]^2, na.rm = TRUE)/2
    return(L)
}

#' @title
#' Gradient supplementary idividuals
#' @param xs Input matrix supplematary individuals.
#' @param B Resamples number.
#' @param k Dimensions number.
#' @param par Parameters.
#' @param lambda Penalization parameter.
#' @param ... other arguments
#' @export
Indsup.GradBIN <- function(xs, B, par, k, lambda, ...) {
    x = as.matrix(xs)
    Bs = as.matrix(B)
    n = nrow(x)
    p = ncol(x)
    aik = n * k
    dTheta = aik

    A = matrix(par[1:aik], n, k)

    lin = cbind(rep(1, n), A) %*% t(Bs)
    hX = exp(lin)/(1 + exp(lin))

    dA = (hX - x) %*% Bs[, -1] + lambda * A
    gradient = c(dA)
    return(gradient)
}


#' @importFrom stats plogis
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
    hX = plogis(lin)

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
    hX = plogis(lin)

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

    hX = plogis(lin)

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
    hX = plogis(lin)

    dA = (hX - x) %*% Bs[, -1] + lambda * A
    gradient = c(dA)
    return(gradient)
}

#' @title
#' missing pattern
#' @param data data frame
#' @param K number of folds

train_miss_pattern <- function(data, K = 7){
    X <- as.matrix(data)
    out <- sapply(1:K, function(x) seq(x, length(X), by = K), simplify = FALSE)
    train <- list()
    for(i in 1:K){
        temp <- X
        temp[out[[i]]] <- NA
        train[[i]] <- temp
    }

    lout <- list(missWold = out, Xtrain = train)
    return(lout)
}


thresholds <- function(x, P, ncuts = 100){
    P <- as.matrix(P)
    if(is.null(colnames(x))) colnames(x) <- paste0('V', 1:ncol(x))

    N1 <- ceiling(sum(x, na.rm = T))
    N0 <- length(as.matrix(x)) - N1
    TE <- lapply(seq(0, 1, length.out = ncuts), function(z){
        Pr <- ifelse(P>=z, 1, 0)
        c1 <- 1 - apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE)/apply(x == 1, 2, sum)
        c2 <- 1 - apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE)/apply(x == 0, 2, sum)
        TE <- 100/2 * (c1 + c2)
        lista <- list(TE)

    })
    TEp <-  data.frame(dplyr::bind_rows(TE), threshold = seq(0, 1, length.out = ncuts))

    thresholds <- TEp %>%
        tidyr::pivot_longer(-threshold, names_to = "variable", values_to = "BACC") %>%
        group_by(variable) %>%
        mutate(merror = min(BACC)) %>%
        dplyr::filter(BACC == merror) %>%
        mutate(row = dplyr::row_number()) %>%
        filter(row == 1) %>% ungroup %>%
        dplyr::select(variable, threshold, BACC)

    thresholds <- thresholds[match(colnames(x), thresholds$variable),]
    Pr <- matrix(NA, nrow(P), ncol(P))
    for(p in 1:ncol(P)){
        Pr[,p] <- ifelse( P[,p] >= thresholds$threshold[p], 1, 0)
    }
    c1 <- 1 - sum(apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE))/N1
    c2 <- 1 - sum(apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE))/N0
    BACC <- round(100/2 * (c1 + c2), 2)

    out <- list(pred = Pr, thres = thresholds, BACC = BACC)
}

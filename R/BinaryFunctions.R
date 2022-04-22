#' @importFrom stats plogis

J.BipLog.BIN <- function(xt, par, k, lambda, ...) {
    #---- Para iniciar.
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

    J = -sum(x * log(hX) + (1 - x) * log(1 - hX), na.rm = TRUE) + lambda *
        sum(A^2, na.rm = TRUE)/2 + lambda * sum(B[, -1]^2, na.rm = TRUE)/2
    return(J)
}

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

    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("Package \"tidyr\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Package \"dplyr\" needed for this function to work. Please install it.",
             call. = FALSE)
    }

    P <- as.matrix(P)
    if(is.null(colnames(x))) colnames(x) <- paste0('V', 1:ncol(x))

    N1 <- ceiling(sum(x, na.rm = T))
    N0 <- length(as.matrix(x)) - N1 - sum(is.na(x))
    TE <- lapply(seq(0, 1, length.out = ncuts), function(z){
        Pr <- ifelse(P>=z, 1, 0)
        c1 <- 1 - apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE)/apply(x == 1, 2, sum, na.rm = T)
        c2 <- 1 - apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE)/apply(x == 0, 2, sum, na.rm = T)
        TE <- 100/2 * (c1 + c2)
        lista <- list(TE)

    })
    TEp <-  data.frame(dplyr::bind_rows(TE), threshold = seq(0, 1, length.out = ncuts))

    thresholds <- TEp |>
        tidyr::pivot_longer(-threshold, names_to = "variable", values_to = "BACC") |>
        dplyr::group_by(variable) |>
        dplyr::mutate(merror = min(BACC)) |>
        dplyr::filter(BACC == merror) |>
        dplyr::mutate(row = dplyr::row_number(), nclass = max(row)) |>
        dplyr::mutate(sel = ifelse(nclass > 1 & row ==2, 1,
                                   ifelse(nclass == 1 & row ==1, 1, 0))) |>
        dplyr::filter(sel == 1) |> dplyr::ungroup() |>
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


crossval <- function(x, k = 2, K = 7, thres = NULL, method = NULL, type = NULL){

    folds <- train_miss_pattern(x, K = K)
    missVal <- folds$missWold
    train <- folds$Xtrain

    err_T = list()
    cv_errD = list()
    for(i in 1:K){
        if(method == "CG" & k > 0){
            missBip <- BiplotML::LogBip(train[[i]], k = k,
                                        method = method, type = type, plot = FALSE)
            P <- fitted_LB(missBip, type = "response")
        }else if(method == "BFGS" & k > 0){
            missBip <- BiplotML::LogBip(train[[i]], k = k,
                                        method = method, plot = FALSE)
            P <- fitted_LB(missBip, type = "response")
        }else if(method == "MM" & k > 0){
            missBip <- BiplotML::LogBip(train[[i]], k = k,
                                        method = "MM", random_start = FALSE, plot = FALSE)
            P <- fitted_LB(missBip, type = "response")
        }else{
            P <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(train[[i]], na.rm=TRUE)))
        }

        Xhat <- matrix(NA, nrow(P), ncol(P))
        for(p in 1:ncol(P)){
            Xhat[,p] <- ifelse(P[,p] >= thres$threshold[p], 1, 0)
        }
        Xhat[is.na(Xhat)] <- 0

        Xhat_pred <- Xhat[missVal[[i]]]
        xReal <- x[missVal[[i]]]

        n1 <- sum(xReal);            n1t <- sum(x)
        n0 <- length(xReal) - n1;    n0t <- length(x) - n1t

        err0 <- sum(ifelse(xReal == 0 & Xhat_pred == 1, 1, 0))/n0
        err1 <- sum(ifelse(xReal == 1 & Xhat_pred == 0, 1, 0))/n1

        err0t <- sum(ifelse(x == 0 & Xhat == 1, 1, 0))/n0t
        err1t <- sum(ifelse(x == 1 & Xhat == 0, 1, 0))/n1t

        err_T[[i]] = 100/2 * (err0t + err1t)

        cv_errD[[i]] <- 100/2 * (err0 + err1)
    }
    cvT <- round(mean(sapply(err_T, mean, na.rm = TRUE), na.rm = TRUE), 2)
    cvD <- round(mean(sapply(cv_errD, mean, na.rm = TRUE), na.rm = TRUE), 2)

    out <- list(cvT = cvT, cvD = cvD)
    return(out)
}

log_like <- function(x, w = NULL, theta){
  n <- nrow(x)
  p <- ncol(x)

  if (!is.null(w)) {
    w = matrix(rep(1, n*p), n, p)
  }
  e_log_like <- matrix(NA, nrow = n, ncol = p)
  P <- plogis(theta)

  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      e_log_like[i,j] <- w[i,j] * (x[i,j] * log(P[i, j]) + (1 - x[i,j]) * log(1 - P[i, j]))
    }
  }
  log_like <- -sum(e_log_like)
  return(log_like)
}

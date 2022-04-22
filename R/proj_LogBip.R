#' @export
#' @title
#' Fitting a Binary Logistic Biplot with Missing Data Using Data Projection and a Block Coordinate Descending Algorithm
#' @description
#' This function impute the missing values of a binary dataset \eqn{X}, and estimates the vector \eqn{\mu}, matrix A and matrix B using data projection model with a block coordinate descending algorithm.
#' @return
#' Imputed \eqn{X} matrix and coordenates of the matrix A and B, and \eqn{\mu}
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x binary matrix.
#' @param k dimensions number. By default \code{k = 2}.
#' @param max_iters maximum iterations.
#' @param random_start random initialization
#' @param epsilon convergence criteria
#' @references
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2022). Logistic biplot with missing data.
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2021). Logistic Biplot by Conjugate Gradient Algorithms and Iterated SVD. Mathematics, 9(16).
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{cv_LogBip}}
#' @examples
#' \donttest{
#' data("Methylation")
#' set.seed(12345)
#' n <- nrow(Methylation)
#' p <- ncol(Methylation)
#' miss <- matrix(rbinom(n*p, 1, 0.2), n, p) #I simulate some missing data
#' miss <- ifelse(miss == 1, NA, miss)
#' x <- Methylation + miss  #Matrix containing missing data
#' out <- LogBip(x, method = "PDLB", maxit = 1000)
#' }

proj_LogBip <- function(x, k = 2, max_iters = 1000,
                    random_start = FALSE, epsilon = 1e-05){

  x <- as.matrix(x)

  verify <- apply(x, 2, sd, na.rm=T)

  if(any(verify == 0)) {
    stop("Some variables have zero variance, so the procedure cannot be applied.")
  }

  n=nrow(x); p=ncol(x);

  partial_decomp = FALSE
  if(n > 100 | p > 20){
    partial_decomp = TRUE
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      message("RSpectra must be installed to optimize the algorithm")
      partial_decomp = FALSE
    }
  }

  if(k > p){
    warning("k must be less than or equal to p. Setting k = p")
    k = p
  }

  W <- ifelse(!is.na(x), 1, 0)
  xi <- ifelse(is.na(x), 1, x)
  x0 <- W * xi + 0.5 * (1 - W)


  if(random_start){
    b0 <- rep(0, p)
    V <- matrix(rnorm(p * k), p, k)
  }else{
    b0 <- colMeans(x0, na.rm = TRUE)
    if (partial_decomp) {
      udv  <- RSpectra::svds(scale(x0, center = b0, scale = FALSE), k = k)
    } else {
      udv = svd(scale(x0, center = b0, scale = FALSE))
    }
    V = matrix(udv$v[, 1:k], p, k)
  }

  A = (x0 - rep(1, n) %*% t(b0)) %*% V; B = V
  theta = rep(1, n) %*% t(b0) + A%*%t(B)

  P = plogis(theta)

  thresh <- thresholds(x0, P, ncuts = 100)$thres

  err <- numeric(max_iters + 1)

  xl <- x0
  err[1] <- log_like(x = xl, w = W, theta = theta)/(n*p)

  for (l in 1:max_iters) {
    Zl <- theta + 4 * (xl - P)
    Ml <- W * Zl + (1 - W) * theta
    Xlc <- xl - rep(1, n)%*%t(b0)
    Mlc <- Ml - rep(1, n)%*%t(b0)
    Yl <- -(t(Xlc) %*% Xlc - t(Xlc) %*% Mlc - t(Mlc) %*% Xlc)
    sdv <- RSpectra::eigs_sym(Yl, k, which = "LM")
    if (any(sdv$values[1:k] < 0)) {
      sdv <- eigen(Yl, symmetric = TRUE)
    }

    V <- matrix(sdv$vectors[, 1:k], p, k)

    b0 <- 1/n * t(Ml - xl %*% tcrossprod(V)) %*% rep(1,n)
    theta <- rep(1, n) %*% t(b0) + (xl - rep(1, n) %*% t(b0)) %*% tcrossprod(V)
    P <- plogis(theta)

    alphas <- rep(1, n) %*% t(as.matrix(thresh$threshold))
    xhat <- ifelse(P < alphas, 0, 1)
    xl <- W * xi + (1 - W) * xhat

    if (l == max_iters) {
      warning("Algorithm has reached ", max_iters, " iterations without converging.")
    }

    err[l + 1] <- log_like(x = xl, w = W, theta = theta) / (dim(x)[1]*dim(x)[2])

    if (l > 2) {
      if ((err[l] - err[l + 1])/err[l] < epsilon)
        break
    }
  }
  A <- (xl - rep(1, n) %*% t(b0)) %*% V
  B <- V

  out <-  list(mu = b0, A = A, B = B, x_est = xl, iter = l, loss_funct = err[1:l])
  return(out)
}




#' @importFrom stats qlogis
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @export
#'
#' @title
#' Multivariate binary data
#' @description
#' Simulate a binary data matrix based on a latent variables model
#' @return
#' X: binary matrix,
#' P: predicted matrix,
#' Theta: matrix of natural parameters,
#' A: row markers,
#' B: column markers,
#' mu: offset term,
#' D: sparsity level,
#' n: number of rows,
#' p: number of columns
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param n number of rows
#' @param p number of columns
#' @param k number of underlying dimensions in the model
#' @param D sparsity control
#' @param C variance control
#' @seealso \code{\link{cv_LogBip}}
#' @examples
#' x <- simBin(n = 100, p = 50, k = 3, D = 0.5)

simBin <- function(n, p, k, D, C = 1){

  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package \"pracma\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package \"mvtnorm\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  Bp <- matrix(rnorm(p*k), p)
  B <- pracma::gramSchmidt(Bp)$Q

  S <-  diag(k)
  mu <- c(rep(qlogis(D), p))

  centro <- c(0, 0, 0)
  A <- mvtnorm::rmvnorm(n, mean = centro, sigma = S)

  logOdds <- rep(1,n)%*%t(mu) + C *(A %*% t(B))
  P <- plogis(logOdds)

  M <- matrix(NA, n, p)
  X <- t(sapply(1:n, function(i) sapply(1:p, function(j) M[i, j] = rbinom(1, 1, P[i,j]))))

  desb <- sum(X)/length(X)

  out <- list(X = X, P = P, Theta = logOdds, A = A, B = B, mu = mu, D = desb, n = n, p = p)

}

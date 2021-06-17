#' @export
#' @title
#' Fitting a Binary Logistic Biplot using coordinate descendent MM algorithm
#' @description
#' This function estimates the vector \eqn{\mu}, matrix A and matrix B using coordinate descendent MM algorithm.
#' @return
#' Coordenates of the matrix A and B, and \eqn{\mu}
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x binary matrix.
#' @param k dimensions number. By default \code{k = 2}.
#' @param iterations maximum iterations.
#' @param truncated if TRUE, find the k largest singular values and vectors of a matrix.
#' @param random random initialization
#' @param epsilon convergence criteria

sdv_MM <- function(x, k = 2, iterations = 1000, truncated = TRUE,
                         random = FALSE, epsilon = 1e-4){
  n <- nrow(x); p <- ncol(x)

  if(k > p){
    message("The value of k must be less than the number of columns in the matrix")
  }

  if(!random){
    mu <- colMeans(8 * x)

    if(truncated){
      vp <- RSpectra::svds(scale(8 * x, center = TRUE, scale = FALSE), k = k, nu = k, nv = k)
    }else{
      vp <- svd(scale(8 * x, center = TRUE, scale = FALSE))
    }

    if(k == 1){
       A <- matrix(vp$u[, 1:k], n, k) * vp$d[1]
       B <- matrix(vp$v[, 1:k], p, k)
    }else{
       A <- vp$u[, 1:k] %*% diag(vp$d[1:k])
       B <- vp$v[, 1:k]
    }
  }else{
    mu <- runif(p)
    A <- matrix(runif(n * k), n, k)
    B <- matrix(runif(p * k), p, k)
  }
  theta <- rep(1,n) %*% t(mu) + (A %*% t(B))

  loss_func <- numeric(iterations + 1)
  pi <- plogis(theta)
  loss_func[1] <- -sum(x * log(pi) + (1 - x) * log(1 - pi), na.rm = TRUE) / (n*p)

  for (i in 1:iterations) {
    old_mu <- mu
    old_A <- A
    old_B <- B

    Z <- theta + 4 * (x - pi)
    mu <- colMeans(Z)

    if(truncated){
      vp <- RSpectra::svds(scale(Z, center = TRUE, scale = FALSE), k = k, nu = k, nv = k)
    }else{
      vp <- svd(scale(Z, center = TRUE, scale = FALSE))
    }
    if(k == 1){
      A <- matrix(vp$u[, 1:k], n, k) * vp$d[1]
      B <- matrix(vp$v[, 1:k], p, k)
    }else{
      A <- vp$u[, 1:k] %*% diag(vp$d[1:k])
      B <- vp$v[, 1:k]
    }
    theta <- rep(1,n) %*% t(mu) + (A %*% t(B))
    pi <- plogis(theta)
    loss_func[i + 1] <- -sum(x * log(pi) + (1 - x) * log(1 - pi), na.rm = TRUE) / (n*p)

    if(i > 10){
      if((loss_func[i] - loss_func[i + 1])/loss_func[i] < epsilon)
        break
    }
    if(i == iterations){
      warning("The algorithm has reached the limit ", iterations, " of iterations without converging")
    }

    if(loss_func[i] < loss_func[i + 1]){
      mu <- old_mu; A <- old_A; B <- old_B; i <- i - 1
      if(truncated){
        warning("The algorithm stopped because the loss function was increased in the iteration ", i, " rerun the algorithm using the argument truncated = FALSE")
      }
    }
    j=i
  }
  out <- list(mu = mu, A = A, B = B, iterations = j, loss_func = loss_func[1:j])
  return(out)
}

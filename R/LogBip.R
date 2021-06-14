#' @importFrom optimr optimr
#' @importFrom stats runif
#' @importFrom logisticPCA logisticSVD
#' @export
#'
#' @title
#' Fitting a Binary Logistic Biplot using bootstrap methodology and select a threshold to create the classification rule
#' @description
#' This function estimates the vector \eqn{\mu}, matrix A and matrix B using the optimization algorithm chosen by the user and applies a bootstrap methodology to determine the confidence ellipses.
#' @return
#' Coordenates of the matrix A and B, threshold for classification rule
#' @details
#' The methods that can be used to estimate the parameters of a logistic biplot
#'
#' - For methods based on the conjugate gradient use method = "CG" and
#'
#'       - type = 1 for the Fletcher Reeves.
#'       - type = 2 for Polak Ribiere.
#'       - type = 3 for Hestenes Stiefel.
#'       - type = 4 for Dai Yuan.
#'
#' - To use the iterative algorithm Majorization-Minimization then method = "MM".
#'
#' - To use the BFGS formula, method = "BFGS".
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x Binary matrix.
#' @param k Dimensions number. By default \code{k = 2}.
#' @param L Penalization parameter. By default \code{L = 0}.
#' @param method Method to be used to estimate the parameters. By default \code{method="CG"}
#' @param type For the conjugate-gradients method. Takes value 1 for the Fletcher–Reeves update, 2 for Polak–Ribiere and 3 for Beale–Sorenson.
#' @param plot Plot the Bootstrap Logistic Biplot.
#' @param maxit The maximum number of iterations. Defaults to 100 for the gradient methods, and 500 without gradient.
#' @param endsegm The segment starts at 0.5 and ends at this value. By default \code{endsegm = 0.90}.
#' @param label.ind By default the row points are not labelled.
#' @param col.ind Color for the rows marks.
#' @param draw The graph to draw ("ind" for the individuals, "var" for the variables and "biplot" for the row and columns coordinates in the same graph)
#' @param random_start Logical value; whether to randomly inititalize the parameters. If \code{FALSE},
#'   algorithm will use an SVD as starting value.
#' @references
#' De Leeuw, Jan (2006). Principal component analysis of binary data by iterated singular value decomposition. Computational Statistics & Data Analysis 50 (1), 21–-39
#'
#' John C. Nash (2011). Unifying Optimization Algorithms to Aid Software System Users:optimx for R. Journal of Statistical Software. 43(9). 1--14.
#'
#' John C. Nash (2014). On Best Practice Optimization Methods in R. Journal of Statistical Software. 60(2). 1--14.
#'
#' Landgraf, A. J., & Lee, Y. (2020). Dimensionality reduction for binary data through the projection of natural parameters. Journal of Multivariate Analysis, 104668.
#'
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{plotBLB}, \link{pred_LB}, \link{fitted_LB}}
#' @examples
#' \dontrun{
#' data("Methylation")
#' ### Conjugate Gradient with Fletcher and Reeves method
#' res <- LogBip(x = Methylation, plot = FALSE)
#' ### Conjugate Gradient with Polak Ribiere method
#' res <- LogBip(x = Methylation, type = 2)
#' ### Majorization-Minimization method
#' res <- LogBip(x = Methylation, method = "MM", maxit = 1000)
#' ### Quasi-Newton algorithm with BFGS method
#' res <- LogBip(x = Methylation, method = "BFGS")
#' }

LogBip <- function(x, k=2, L=0, method="CG", type = 1, plot=TRUE, maxit=NULL, endsegm = 0.90, label.ind = FALSE, col.ind = NULL,
                   draw = c("biplot","ind","var"), random_start=FALSE){

  if(any(is.na(x))){
    x <-  sweep(x, MARGIN = 2,
                 STATS = colMeans(x, na.rm = TRUE),
                 FUN =  function(z,s) ifelse(is.na(z), s, z)
    )}

  n=nrow(x); p=ncol(x); aik=n*k; bjk=p*(k+1)
  dTheta = aik + bjk; s=k+1

  partial_decomp = TRUE
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    message("RSpectra must be installed to optimize the algorithm")
    partial_decomp = FALSE
  }

  if(random_start){
    params <- runif(dTheta)
  }else{
    b0 <- colMeans(4*x, na.rm = TRUE)
    if (partial_decomp) {
      udv  <- RSpectra::svds(scale(4*x, center = b0, scale = FALSE), k = k)
    } else {
      udv = svd(scale(4*x, center = b0, scale = FALSE))
    }
    A <- matrix(udv$u[, 1:k], n, k) %*% diag(udv$d[1:k], nrow = k, ncol = k)
    B <- matrix(udv$v[, 1:k], p, k)

    params <- c(A, b0, B)
  }

  if(method == "CG"){
    res <- optimr(par=params, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
                 xt=x, k = k, lambda = L, method = method, control = list(type = type))
  }else if(method == "MM"){
    if(!is.null(maxit)){
    res <- logisticPCA::logisticSVD(x = x, k = k, max_iters = maxit, random_start = random_start)
    }else{
    res <- logisticPCA::logisticSVD(x = x, k = k, max_iters = 2000, random_start = random_start)
    }
  }else{
    res <- optimr(par=params, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
                 xt=x, k = k, lambda = L, method = method)
  }

  ### row coordenates
  if(method == "MM"){
    Ahat <- data.frame(res$A)
  }else{
  par <- res$par
  Ahat <- data.frame(matrix(res$par[1:aik], n, k))
  }
  colnames(Ahat) = c(paste0("Dim", seq(1,k,1)))
  rownames(Ahat) = rownames(x)

  ### column coordinates
  if(method == "MM"){
  Bhat <- data.frame(res$mu, res$B)
  }else{
  Bhat <- data.frame(matrix(par[(aik + 1):dTheta], p, k+1))
  }
  rownames(Bhat) = colnames(x)
  colnames(Bhat) = c(paste0("bb", seq(0,k,1)))

  rownames(Ahat) <- rownames(x)

  if(method == "MM") method <- "MM: Majorization-Minimization"
  if(method == "CG" & type == 1) method <- "CG: Fletcher--Reeves"
  if(method == "CG" & type == 2) method <- "CG: Polak--Ribiere"
  if(method == "CG" & type == 3) method <- "CG: Beale--Sorenson"
  if(method == "Rcgmin") method <- "CG: Dai--Yuan"

  out <- list(Ahat = Ahat, Bhat = Bhat, method=method)

  if (plot & ncol(Ahat)>1) {
    print(plotBLB(x=out, ellipses = FALSE, endsegm = endsegm, titles = "Logistic Biplot", label.ind = label.ind, col.ind = col.ind, draw = draw))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

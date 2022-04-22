#' @importFrom optimr optimr
#' @importFrom stats runif
#' @importFrom stats sd
#' @export
#' @title
#' Fitting a Binary Logistic Biplot using optimization methods
#' @description
#' This function estimates the vector \eqn{\mu}, matrix A and matrix B using the optimization algorithm chosen by the user.
#' The PDLB method allows to enter a binary matrix with missing data
#' @return
#' Coordenates of the matrix A and B, threshold for classification rule.
#' Furthemore, for the PDLB method, the imputed matrix is returned.
#' @details
#' The methods that can be used to estimate the parameters of a logistic biplot
#'
#' - For methods based on the conjugate gradient use method = "CG" and
#'
#'       type = 1 for the Fletcher Reeves; type = 2 for Polak Ribiere; type = 3 for Hestenes Stiefel and type = 4 for Dai Yuan.
#'
#' - To use the iterative coordinate descendent MM algorithm then method = "MM".
#'
#' - If the binary matrix X has missing data, use method = "PDLB". In case it's required to estimate the row coordinates of other individuals, this method is also the most appropriate. For more details see the paper "Logistic biplot with missing data".
#'
#' - To use the BFGS formula, method = "BFGS".
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x Binary matrix.
#' @param k Dimensions number. By default \code{k = 2}.
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
#' @param L Penalization parameter. By default \code{L = 0}.
#' @param cv_LogBip Indicates if the procedure is being used for cross validation.
#' @references
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2022). Logistic biplot with missing data. In Process.
#'
#' Babativa-Marquez, J. G., & Vicente-Villardon, J. L. (2021). Logistic Biplot by Conjugate Gradient Algorithms and Iterated SVD. Mathematics, 9(16).
#'
#' John C. Nash (2011). Unifying Optimization Algorithms to Aid Software System Users:optimx for R. Journal of Statistical Software. 43(9). 1--14.
#'
#' John C. Nash (2014). On Best Practice Optimization Methods in R. Journal of Statistical Software. 60(2). 1--14.
#'
#' Nocedal, J.;Wright, S. (2006). Numerical optimization; Springer Science & Business Media.
#'
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{plotBLB}, \link{pred_LB}, \link{fitted_LB}}
#' @examples
#' \donttest{
#' data("Methylation")
#' # If the binary matrix has no missing data and does not require the projection
#' # of supplementary individuals, you can use an coordinate descendent MM algorithm
#' res_MM <- LogBip(x = Methylation, method = "MM", maxit = 1000)
#' # If the binary matrix has missing data or requires the projection of supplementary
#' #individuals, use a method based on data projection with a block coordinate descent algorithm
#' data("Methylation")
#' set.seed(12345)
#' n <- nrow(Methylation)
#' p <- ncol(Methylation)
#' miss <- matrix(rbinom(n*p, 1, 0.2), n, p) #I simulate some missing data
#' miss <- ifelse(miss == 1, NA, miss)
#' x <- Methylation + miss  #Matrix containing missing data
#' out <- LogBip(x, method = "PDLB", maxit = 1000)
#' }

LogBip <- function(x, k=2, method="MM", type = NULL, plot=TRUE, maxit=NULL, endsegm = 0.90, label.ind = FALSE, col.ind = NULL,
                   draw = c("biplot","ind","var"), random_start=FALSE, L = 0, cv_LogBip = FALSE){

  x <- as.matrix(x)
  verify <- apply(x, 2, sd, na.rm=T)

  if(any(verify == 0)) {
    stop("Some variables have zero variance, so the procedure cannot be applied.")
  }


  if(cv_LogBip & any(is.na(x))){
    x <-  sweep(x, MARGIN = 2,
                STATS = colMeans(x, na.rm = TRUE),
                FUN =  function(z,s) ifelse(is.na(z), s, z)
    )}else if(any(is.na(x)) & method != "PDLB"){
    warning("Binary matrix contains missing values, the ", method, " method has been changed to PDLB method.")
    method = "PDLB"
  }

  n=nrow(x); p=ncol(x); aik=n*k; bjk=p*(k+1)
  dTheta = aik + bjk; s=k+1

  truncated = FALSE
  if(n > 100 | p > 20){
    truncated = TRUE
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      message("RSpectra must be installed to optimize the algorithm")
      truncated = FALSE
    }
  }

  if(method != "PDLB"){
    if(random_start){
      params <- runif(dTheta)
    }else{
      b0 <- colMeans(4*x, na.rm = TRUE)
      if (truncated) {
        udv  <- RSpectra::svds(scale(4*x, center = b0, scale = FALSE), k = k)
      } else {
        udv = svd(scale(4*x, center = b0, scale = FALSE))
      }
      A <- matrix(udv$u[, 1:k], n, k) %*% diag(udv$d[1:k], nrow = k, ncol = k)
      B <- matrix(udv$v[, 1:k], p, k)

      params <- c(A, b0, B)
    }
  }

  if(method == "CG"){
    res <- optimr(par=params, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
                 xt=x, k = k, lambda = L, method = method, control = list(type = type))
  }else if(method == "MM"){
    if(!is.null(maxit)){
    res <- sdv_MM(x = x, k = k, iterations = maxit, random = random_start, truncated = truncated)
    }else{
    res <- sdv_MM(x = x, k = k, random = random_start, truncated = truncated)
    }
  }else if(method == "PDLB"){
    if(!is.null(maxit)){
      res <- proj_LogBip(x = x, k = k, max_iters = maxit, random_start = random_start)
    }else{
      res <- proj_LogBip(x = x, k = k, random_start = random_start)
    }
  }else{
    res <- optimr(par=params, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
                 xt=x, k = k, lambda = L, method = method)
  }

  if(method %in% c("MM", "PDLB")){
    Ahat <- data.frame(res$A)
  }else{
  par <- res$par
  Ahat <- data.frame(matrix(res$par[1:aik], n, k))
  }
  colnames(Ahat) = c(paste0("Dim", seq(1,k,1)))
  rownames(Ahat) = rownames(x)

  if(method %in% c("MM", "PDLB")){
  Bhat <- data.frame(res$mu, res$B)
  }else{
  Bhat <- data.frame(matrix(par[(aik + 1):dTheta], p, k+1))
  }
  rownames(Bhat) = colnames(x)
  colnames(Bhat) = c(paste0("bb", seq(0,k,1)))

  if(method == "MM"){
    method <- "coordinate descendent MM"
    out <- list(Ahat = Ahat, Bhat = Bhat, method=method,
                loss_function = res$loss_func, iterations = res$iterations)
  }else if(method == "PDLB"){
    method <- "data projection using a block coordinate descent"
    out <- list(Ahat = Ahat, Bhat = Bhat, impute_x = res$x_est, method=method,
                loss_function = res$loss_funct, iterations = res$iter)
  }else{
  if(method == "CG" & type == 1) method <- "CG: Fletcher--Reeves"
  if(method == "CG" & type == 2) method <- "CG: Polak--Ribiere"
  if(method == "CG" & type == 3) method <- "CG: Beale--Sorenson"
  if(method == "Rcgmin") method <- "CG: Dai--Yuan"
  out <- list(Ahat = Ahat, Bhat = Bhat, method=method)
  }


  if (plot & ncol(Ahat)>1) {
    print(plotBLB(x=out, ellipses = FALSE, endsegm = endsegm, titles = "Logistic Biplot", label.ind = label.ind, col.ind = col.ind, draw = draw))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

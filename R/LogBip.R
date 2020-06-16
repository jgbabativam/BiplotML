#' @importFrom optimr optimr
#' @importFrom stats runif
#' @export
#'
#' @title
#' Fitting a Binary Logistic Biplot using bootstrap methodology
#' @description
#' This function estimates the vector \eqn{\mu}, matrix A and matrix B using the optimization algorithm chosen by the user and applies a bootstrap methodology to determine the confidence ellipses.
#' @return
#' Coordenates of the matrix A and B in resamples and Biplot
#' @details
#' Fitting when sup=TRUE ... whereas sup=FALSE ...
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x Binary matrix.
#' @param k Dimensions number. By default \code{k = 2}.
#' @param L Penalization parameter. By default \code{L = 0}.
#' @param method Method to be used to estimate the parameters. By default \code{ method="CG"}
#' @param type For the conjugate-gradients method. Takes value 1 for the Fletcher–Reeves update, 2 for Polak–Ribiere and 3 for Beale–Sorenson.
#' @param plot Plot the Bootstrap Logistic Biplot.
#' @param maxit The maximum number of iterations. Defaults to 100 for the gradient methods, and 500 without gradient.
#' @param endsegm The segment starts at 0.5 and ends at this value. By default \code{endsegm = 0.90}.
#' @param label.ind By default the row points are not labelled.
#' @param col.ind Color for the rows marks.
#' @param draw The graph to draw ("ind" for the individuals, "var" for the variables and "biplot" for the row and columns coordinates in the same graph)
#' @references
#' John C. Nash (2011). Unifying Optimization Algorithms to Aid Software System Users:optimx for R. Journal of Statistical Software. 43(9). 1--14.
#'
#' John C. Nash (2014). On Best Practice Optimization Methods in R. Journal of Statistical Software. 60(2). 1--14.
#'
#' Milan, L., & Whittaker, J. (1995). Application of the parametric bootstrap to models that incorporate a singular value decomposition. Applied Statistics, 44, 31–49.
#'
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{plotBLB}, \link{performanceBLB}}
#' @examples
#' data("Methylation")
#' set.seed(02052020)
#' ### Conjugate Gradient with Fletcher and Reeves method
#' res.LB <- LogBip(x = Methylation)
#' ### Quasi-Newton algorithm with BFGS method
#' res.LB <- LogBip(x = Methylation, method = "BFGS")

LogBip <- function(x, k=2, L=0, method="CG", type = 1, plot=TRUE, maxit=NULL, endsegm = 0.90, label.ind = FALSE, col.ind = NULL,
                   draw = c("biplot","ind","var")){
  n=nrow(x); p=ncol(x); aik=n*k; bjk=p*(k+1)
  dTheta = aik + bjk; s=k+1

  if(method == "CG"){
    res = optimr(par=runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
                 xt=x, k = k, lambda = L, method = method, control = list(type = type))
  }else{
    res = optimr(par=runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
                 xt=x, k = k, lambda = L, method = method)
  }

  par=res$par
  ### row coordenates
  Ahat = data.frame(matrix(res$par[1:aik], n, k))
  colnames(Ahat) = c(paste0("Dim", seq(1,k,1)))
  rownames(Ahat) = rownames(x)

  ### column coordinates
  Bhat = data.frame(matrix(par[(aik + 1):dTheta], p, k+1))
  rownames(Bhat) = colnames(x)
  colnames(Bhat) = c(paste0("bb", seq(0,k,1)))

  LogP = as.matrix(cbind(rep(1,nrow(x)),Ahat))%*%t(Bhat)
  P = exp(LogP)/(1+exp(LogP))
  Pr = ifelse(P>=0.5, 1, 0)

  PCC = ifelse((x==1 & Pr==1) | (x==0 & Pr==0), 1, 0)
  PCC[is.na(PCC)] <- 0
  ones <-  apply(x, 2, sum, na.rm=TRUE)
  zeros <- n - ones

  confusion <- data.frame( Sensitivy = round(100*apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE)/ ones, 1),
                           Specificity = round(100*apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE)/ zeros, 1),
                           Global = round(100*colSums(PCC)/nrow(PCC), 1))

  rownames(Ahat) <- rownames(x)
  rownames(Bhat) <- colnames(x)

  if(method == "CG" & type == 1) method <- "CG: Fletcher--Reeves"
  if(method == "CG" & type == 2) method <- "CG: Polak--Ribiere"
  if(method == "CG" & type == 3) method <- "CG: Beale--Sorenson"
  if(method == "Rcgmin") method <- "CG: Dai--Yuan"

  out <- list(Ahat = Ahat, Bhat = Bhat, pred= Pr, fit = confusion, method=method)

  if (plot & ncol(Ahat)>1) {
    print(plotBLB(x=out, ellipses = FALSE, endsegm = endsegm, titles = "Logistic Biplot", label.ind = label.ind, col.ind = col.ind, draw = draw))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

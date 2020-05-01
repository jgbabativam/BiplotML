#' @importFrom optimr optimr
#' @importFrom shapes procOPA
#' @import dplyr
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
#' @param plot Plot the Bootstrap Logistic Biplot.
#' @param sup Boolean, if TRUE, rows that are not selected in each resample are treated as supplementary individuals. See details.
#' @param ellipses Draw confidence ellipses. By default is TRUE.
#' @param maxit The maximum number of iterations. Defaults to 100 for the gradient methods, and 500 without gradient.
#' @param resamples Number of iterations in the bootstrap process. By default \code{100}.
#' @param conf Level confidence in the ellipses. By default \code{conf=0.90}
#' @references
#' John C. Nash (2011). Unifying Optimization Algorithms to Aid Software System Users:optimx for R. Journal of Statistical Software. 43(9). 1--14.
#'
#' John C. Nash (2014). On Best Practice Optimization Methods in R. Journal of Statistical Software. 60(2). 1--14.
#'
#' Milan, L., & Whittaker, J. (1995). Application of the parametric bootstrap to models that incorporate a singular value decomposition. Applied Statistics, 44, 31–49.
#'
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{plot.BLB}, \link{performanceBLB}}
#' @examples
#' data("Methylation")
#' set.seed(123456)
#' out.sup <- bootBLB(x = Methylation, ellipses = TRUE)
#' out <- bootBLB(x = Methylation, sup = FALSE, ellipses = TRUE)

bootBLB <- function(x, k=2, L=0, method="CG", plot=TRUE, sup=TRUE,
                    ellipses=TRUE, maxit=NULL, resamples = 100, conf = 0.9){
  n=nrow(x); p=ncol(x); aik=n*k; bjk=p*(k+1)
  dTheta = aik + bjk; s=k+1

  res = optimr(par=runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
               xt=x, k = k, lambda = L, method = method)
  par=res$par
  ### row coordenates
  A = data.frame(matrix(res$par[1:aik], n, k))
  colnames(A) = c(paste0("Dim", seq(1,k,1)))
  rownames(A) = rownames(x)

  ### column coordinates
  B = data.frame(matrix(par[(aik + 1):dTheta], p, k+1))
  rownames(B) = colnames(x)
  colnames(B) = c(paste0("b", seq(0,k,1)))

  ### bootstrap process
  Ai = A %>% mutate(rowId = row_number())
  Bi = as.matrix(B[,1:s])

  listA = list()
  listB = list()

  indica = rep(1:nrow(x))
  xtemp <- x %>% `rownames<-`(seq_len(nrow(x)))

  for (i in 1:resamples) {
    sample_ind = as.matrix(sample(indica,replace = T))
    xb = xtemp[sample_ind,]
    n=nrow(xb); p=ncol(xb);
    aik=n*k; bjk=p*(k+1)

    ran1 = as.matrix(runif(aik))
    sample_ind2=rbind(sample_ind,sample_ind+n)

    par1 = ran1[sample_ind2,]
    par2 = runif(bjk)
    param = c(par1, par2)

    res.b = optimr(par=param, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
                   xt=xb, k = k, lambda = L, method = method)

    ### Column Coordinates
    Bb = matrix(res.b$par[(aik + 1):dTheta], p, k+1)

    ### Row coordinates
    Ab = matrix(res.b$par[1:aik], n, k)
    Ab_j = data.frame(Ab, rowId=sample_ind)
    Ab_j = dplyr::distinct(Ab_j, rowId, .keep_all=T)

    if(sup){
      Asup = dplyr::anti_join(Ai, Ab_j, by="rowId")
      n_sup = nrow(Asup)
      pars = runif(dim(Asup)[1] * (dim(Asup)[2]-1))
      xsup = xtemp[Asup$rowId,]

      res.sup = optimr(par=pars, fn = Indsup.BIN, gr = Indsup.GradBIN,
                       xs=xsup, B=Bb, k = k, lambda = L, method = method)

      Asup = data.frame(matrix(res.sup$par, n_sup, k), rowId = Asup[,ncol(Asup)])
      Abr = rbind(Ab_j, Asup) %>%
        arrange(rowId)

      Ab = Abr %>% dplyr::select(-rowId)

      #--- Rotate - reflect matrix A.
      outA = procOPA(as.matrix(A), as.matrix(Ab), scale = TRUE, reflect = TRUE)
      Arot = outA$Bhat
      R = outA$R
      SCR_A = outA$OSS
    }else{
      Abr = Ab_j %>% arrange(rowId)                ### Matriz estimada bootstrap
      Ai_b = inner_join(Ai, Abr, by="rowId") %>%   ### Matriz inicial solo con ind del boot
        dplyr::select(paste0("Dim", seq(1,k,1)))
      Ab = Abr %>% dplyr::select(-rowId)

      #--- Rotate - reflect matrix A.
      outA = procOPA(as.matrix(Ai_b), as.matrix(Ab), scale = TRUE, reflect = TRUE)
      Arot = outA$Bhat
      R = outA$R
      SCR_A = outA$OSS
    }
    Brot = Bb[,2:s] %*% R

    Arot = data.frame(Arot, ind=Abr$rowId)
    colnames(Arot)=c(paste0("Dim", seq(1,k,1)), "ind")
    rownames(Arot)=rownames(x[Abr$rowId,])

    Brot = data.frame(Bb[,1], Brot, param = seq(1,p,1))
    colnames(Brot)=c(paste0("b", seq(0,k,1)), "param")
    rownames(Brot)=rownames(B)

    Arot$resample = i; Brot$resample = i

    listA[[i]] <- Arot
    listB[[i]] <- Brot

  }

  ResultA <- dplyr::bind_rows(listA)
  ResultB <- dplyr::bind_rows(listB)
  rows <- rownames(x)
  cols <- colnames(x)

  # i) Expected values by rows and parameters
  EspA <- aggregate(.~ind, ResultA, mean)
  colnames(EspA) <- c("ind", paste0("Dimb", seq(1,k,1)), "resampleb")

  EspB <- aggregate(.~param, ResultB, mean)
  colnames(EspB) <- c("param", paste0("bb", seq(0,k,1)), "resampleb")

  Ahat <- EspA %>% dplyr::select(starts_with("Dim"))
  Bhat <- EspB %>% dplyr::select(starts_with("bb"))


  if(ellipses){
    # ii.) Center Bootstrap Cloud

    CentBootA <-  left_join(ResultA, EspA, by="ind")%>%
      mutate(Dim1c = Dim1-Dimb1, Dim2c=Dim2-Dimb2) %>%
      dplyr::select(ind, resample, Dim1c, Dim2c)

    CentBootB <-  left_join(ResultB, EspB, by="param") %>%
      mutate(bj0c = b0-bb0, bj1c = b1-bb1, bj2c=b2-bb2) %>%
      dplyr::select(param, resample, bj0c, bj1c, bj2c)

    #...... iii. iv. v.) Se hace para cada individuo
    n = nrow(EspA); p = nrow(EspB)

    datalist = list()
    q <- conf
    for (a in 1:n){
      # El numero de filas no es igual en todos los casos, depende del numero de veces que salio en el
      # proceso de remuestreo
      temp <-  dplyr::filter(CentBootA, ind == a) %>%
        dplyr::select(Dim1c, Dim2c)
      sol <- svd(temp)
      #rep = sol$u %*% diag(sol$d) %*% t(sol$v)
      k <-  sqrt(rowSums(sol$u * sol$u))            # iii.) Distancia de cada punto boots en la base U al centroide
      r <-  quantile(k, q)                          #  iv.) Percentil q%
      theta  <-  seq(0, 2 * pi, length = 500)       #   v.) Angulos para aproximar elipse
      z <-  rbind(r * cos(theta), r * sin(theta))

      #.... Transformación Bootstrap al espacio original
      esp  <-  dplyr::filter(EspA, ind == a) %>%
        dplyr::select(-ind, -resampleb)

      v <-  as.data.frame(as.matrix(rep(1, 500))%*% as.matrix(esp) + r*t(z)%*%diag(sol$d)%*%t(sol$v))
      v$ind  <-  a
      datalist[[a]]  <-  v
    }

    Ellip <- dplyr::bind_rows(datalist)
  }


  LogP = as.matrix(cbind(rep(1,nrow(x)),Ahat))%*%t(Bhat) # Matriz reproducida del Log(P).
  P = exp(LogP)/(1+exp(LogP))                            # Matriz de de probabilidades esperadas.
  Pr = ifelse(P>=0.5, 1, 0)                              # Matriz de valores predichos

  ### Porcentaje de clasificaciones correctas
  PCC = ifelse((x==1 & Pr==1) | (x==0 & Pr==0), 1, 0)
  ones <-  apply(x, 2, sum)
  zeros <- n - ones

  confusion <- data.frame( Sensitivy = round(100*apply((Pr == 1) & (x == 1), 2, sum)/ ones, 1),
                           Specificity = round(100*apply((Pr == 1) & (x == 0), 2, sum)/ zeros, 1),
                           Global = round(100*colSums(PCC)/nrow(PCC), 1))


  rownames(Ahat) <- rownames(x)
  rownames(Bhat) <- colnames(x)

  if(ellipses){
    out <- list(Ahat = Ahat, Bhat = Bhat, BootA=ResultA, BootB=ResultB, pred= Pr, fit = confusion, rows=rows, cols = cols, method=method, Ellip=Ellip)
  }else{
    out <- list(Ahat = Ahat, Bhat = Bhat, BootA=ResultA, BootB=ResultB, pred= Pr, fit = confusion, rows=rows, cols = cols, method=method)
  }

  if (plot & ncol(Ahat)>1) {
    print(plot.BLB(x=out, ellipses = ellipses, endsegm = 0.95))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

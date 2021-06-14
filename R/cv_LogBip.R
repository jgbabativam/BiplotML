#' @export
#' @title
#' Cross-Validation for logistic biplot
#' @description
#' This function run cross validation for logistic biplot
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param data Binary matrix.
#' @param k Dimensions to analyze. By default \code{k = 1:3}.
#' @param K folds. By default \code{K = 7}.
#' @param method Method to be used to estimate the parameters. By default \code{ method="MM"}
#' @param type For the conjugate-gradients method. Takes value 1 for the Fletcher–Reeves update, 2 for Polak–Ribiere and 3 for Beale–Sorenson.
#' @param plot draw the graph. By default \code{plot=TRUE}
#' @param maxit The maximum number of iterations. Defaults to 100 for the gradient methods, and 2000 for MM algorithm.
#'
#' @references
#' Bro R and Kjeldahl K and Smilde AK. (2008). Cross-validation of component models: a critical look at current methods. Analytical and bioanalytical chemistry. 390(5):1241-1251
#'
#' Wold S. (1978). Cross-validatory estimation of the number of components in factor and principal components models. Technometrics. 20(4):397--405.
#'
#' @seealso \code{\link{plotBLB}, \link{pred_LB}, \link{fitted_LB}}
#' @examples
#' data("Methylation")
#' cv_CG <- cv_LogBip(data = Methylation, k = 0:3, method = "CG", type = 1)
#' \dontrun{
#' cv_MM <- cv_LogBip(data = Methylation, k = 0:4, method = "MM", maxit = 1000)
#' }


cv_LogBip <- function(data, k = 0:5, K = 7, method = "MM", type = NULL, plot = TRUE, maxit = NULL){

  x <- as.matrix(data)
  min <- min(k)
  j <- k

  cvD <- matrix(NA, length(k), 3)
  for(k in j){

    if(method == "CG" & k > 0){
      bip <- BiplotML::LogBip(x, k = k, method = method, type = type, plot = FALSE)
      thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds
    }else if(method == "BFGS" & k > 0){
      bip <- BiplotML::LogBip(x, k = k, method = method, plot = FALSE)
      thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds
    }else if(method == "MM" & k > 0){
      if(!is.null(maxit)){
        bip <- BiplotML::LogBip(x, k = k, method = "MM", maxit = maxit, plot = FALSE)
      }else{
        bip <- BiplotML::LogBip(x, k = k, method = "MM", plot = FALSE)
      }
      thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds
    }else{
      P <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(x, na.rm=TRUE)))
      thres <- thresholds(x = x, P = P, ncuts = 50)$thres
    }

    folds <- train_miss_pattern(x, K = K)
    missVal <- folds$missWold
    train <- folds$Xtrain

    cv_errD = list()
    train_err = list()
    for(i in 1:K){
      if(method == "CG" & k > 0){
        missBip <- BiplotML::LogBip(train[[i]], k = k,
                                    method = method, type = type, plot = FALSE)
        P <- fitted_LB(missBip, type = "response")
      }else if(method %in% c("BFGS", "MM") & k > 0){
        if(!is.null(maxit)){
          missBip <- BiplotML::LogBip(train[[i]], k = k,
                                      method = method, maxit = maxit, plot = FALSE)
        }else{
          missBip <- BiplotML::LogBip(train[[i]], k = k,
                                      method = method, plot = FALSE)
        }
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

      n1 <- sum(xReal);          n1t <- sum(x)
      n0 <- length(xReal) - n1;  n0t <- length(x) - n1t

      err0 <- sum(ifelse(xReal == 0 & Xhat_pred == 1, 1, 0))/n0
      err1 <- sum(ifelse(xReal == 1 & Xhat_pred == 0, 1, 0))/n1

      err0t <- sum(ifelse(x == 0 & Xhat == 1, 1, 0))/n0t
      err1t <- sum(ifelse(x == 1 & Xhat == 0, 1, 0))/n1t

      cv_errD[[i]] <- 100/2 * (err0 + err1)
      train_err[[i]] <- 100/2 * (err0t + err1t)
    }
    row <- k - min + 1
    cvD[row, 1] <- k
    cvD[row, 2] <- round(mean(sapply(cv_errD, mean, na.rm = TRUE), na.rm = TRUE), 2)
    cvD[row, 3] <- round(mean(sapply(train_err, mean, na.rm = TRUE), na.rm = TRUE), 2)
  }
  out <- as.data.frame(cvD)
  colnames(out) <- c("k", "cv-error", "train-error")


  if(plot){
    cvp <-   out %>%
      pivot_longer(-k, names_to = "type.error", values_to = "error")   %>%
      ggplot(aes(x = k, y = error, linetype = type.error, shape = type.error, color = type.error)) +
      geom_line() +
      geom_point(size = 3.0) +
      geom_vline(xintercept = out[which.min(out$`cv-error`), "k"], linetype = 2) +
      scale_color_manual(values = c("red", "blue")) +
      labs(x = "dimensions (k)", y = "Classification error (%)")+
      scale_x_continuous(breaks = seq(0,20,1)) +
      theme_classic() + theme(legend.position = "top", legend.title = element_blank(),
                              legend.text=element_text(size=11))
    print(cvp)
  }

  return(out)
}

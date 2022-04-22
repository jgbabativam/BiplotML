#' @export
#' @title
#' Cross-Validation for logistic biplot
#' @description
#' This function run cross-validation for logistic biplot
#' @return
#' Training error and generalization error for a logistic biplot model.
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
#' @seealso \code{\link{LogBip}, \link{pred_LB}, \link{fitted_LB}, \link{simBin}}
#' @examples
#' \donttest{
#' set.seed(1234)
#' x <- simBin(n = 100, p = 50, k = 3, D = 0.5, C = 20)
#' # cross-validation with coordinate descendent MM algorithm
#' cv_MM <- cv_LogBip(data = x$X, k=0:5, method = "MM", maxit = 1000)
#'
#' # cross-validation with CG Fletcher-Reeves algorithm
#' cv_CG <- cv_LogBip(data = x$X, k=0:5, method = "CG", type = 1)
#'
#' # cross-validation with projection data and block coordinate descending algorithm
#' cv_PB <- cv_LogBip(data = x$X, k=0:5, method = "PDLB", maxit = 1000)
#' }

cv_LogBip <- function(data, k = 0:5, K = 7, method = "MM", type = NULL,
                      plot = TRUE, maxit = NULL){

  x <- as.matrix(data)

  verify <- apply(x, 2, sd, na.rm=T)
  if(any(verify == 0)) {
    stop("Some variables have zero variance, so the procedure cannot be applied.")
  }

  min <- min(k)
  j <- k

  if(any(is.na(x)) & method != "PDLB"){
    warning("Binary matrix contains missing values, the ", method, " method has been changed to PDLB method.")
    method = "PDLB"
  }

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
    }else if(method == "PDLB" & k > 0){
      if(!is.null(maxit)){
        bip <- BiplotML::LogBip(x, k = k, method = "PDLB", maxit = maxit, plot = FALSE)
        xt <- bip$impute_x
      }else{
        bip <- BiplotML::LogBip(x, k = k, method = "PDLB", plot = FALSE)
        xt <- bip$impute_x
      }
      thres <- BiplotML::pred_LB(bip, bip$impute_x, ncuts = 50)$thresholds
    }else{
      theta <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(x, na.rm=TRUE)))
      P <- plogis(theta)
      thres <- thresholds(x = x, P = P, ncuts = 50)$thres
    }

    if(k == 0 & any(is.na(x))){
      m <- apply(x, 2, mean, na.rm = T)
      xt <- ifelse(is.na(x) & m > 0.5, 1, 0)
    }

    if(method != "PDLB") xt <- x
    if(k == 0 && any(!is.na(x)) && method == "PDLB") xt <- x

    folds <- train_miss_pattern(xt, K = K)
    missVal <- folds$missWold
    train <- folds$Xtrain

    cv_errD = list()
    train_err = list()
    for(i in 1:K){
      if(method == "CG" & k > 0){
        missBip <- BiplotML::LogBip(train[[i]], k = k,
                                    method = method, type = type, plot = FALSE, cv_LogBip = TRUE)
        P <- fitted_LB(missBip, type = "response")
      }else if(method %in% c("BFGS", "MM", "PDLB") & k > 0){
        if(!is.null(maxit)){
          missBip <- BiplotML::LogBip(train[[i]], k = k,
                                      method = method, maxit = maxit, plot = FALSE, cv_LogBip = TRUE)
        }else{
          missBip <- BiplotML::LogBip(train[[i]], k = k,
                                      method = method, plot = FALSE, cv_LogBip = TRUE)
        }
        P <- fitted_LB(missBip, type = "response")
      }else{
        theta <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(train[[i]], na.rm=TRUE)))
        P <- plogis(theta)
      }

      Xhat <- matrix(NA, nrow(P), ncol(P))
      for(p in 1:ncol(P)){
        Xhat[,p] <- ifelse(P[,p] >= thres$threshold[p], 1, 0)
      }
      Xhat[is.na(Xhat)] <- 0

      Xhat_pred <- Xhat[missVal[[i]]]
      xReal <- xt[missVal[[i]]]

      n1 <- sum(xReal);          n1t <- sum(xt)
      n0 <- length(xReal) - n1;  n0t <- length(xt) - n1t

      err0 <- sum(ifelse(xReal == 0 & Xhat_pred == 1, 1, 0))/n0
      err1 <- sum(ifelse(xReal == 1 & Xhat_pred == 0, 1, 0))/n1

      err0t <- sum(ifelse(xt == 0 & Xhat == 1, 1, 0))/n0t
      err1t <- sum(ifelse(xt == 1 & Xhat == 0, 1, 0))/n1t

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
    cvp <-   out |>
      tidyr::pivot_longer(-k, names_to = "type.error", values_to = "error")   |>
      ggplot2::ggplot(ggplot2::aes(x = k, y = error, linetype = type.error, shape = type.error, color = type.error)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 3.0) +
      ggplot2::geom_vline(xintercept = out[which.min(out$`cv-error`), "k"], linetype = 2) +
      ggplot2::scale_color_manual(values = c("red", "blue")) +
      ggplot2::labs(x = "Dimensions (k)", y = "Classification error (%)",
                    caption = paste("Cross validation using the ", method, " method"))+
      ggplot2::scale_x_continuous(breaks = seq(0,20,1)) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank(),
                              legend.text=ggplot2::element_text(size=11))
    print(cvp)
  }

  return(out)
}

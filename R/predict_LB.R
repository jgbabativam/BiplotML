#' @title
#' Predict logistic biplot and thresholds by variable
#' @description
#' Predicts the binary matrix and calculates the optimal thresholds per variable that minimize the Balanced Accuracy (BACC)
#' @param object BiplotML object
#' @param x Binary matrix.
#' @param ncuts Number of equidistant cuts between 0 and 1 that will be evaluated. By default \code{ncuts = 100}
#' @details
#' The threshold for each variable is lowered to minimize the Balanced Error Rate (BER).
#' \deqn{BACC = \frac{1}{2} (\frac{TP}{TP+FN} + \frac{TN}{TN+FP}),}
#' where \code{TP} is the number of true positives, \code{TN} is the number of true negatives, \code{FP} is the number of false positives and \code{FN} is the number of false negatives
#' @return
#' This function returns the thresholds per variable, the predicted matrix, the confusion matrix and the BER.
#' @examples
#' \donttest{
#' data("Methylation")
#' LB <- LogBip(Methylation, plot = FALSE)
#' out <- pred_LB(LB, Methylation)
#' }
#' @export

pred_LB <- function(object, x, ncuts = 100){

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package \"tidyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package \"dplyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  P <- fitted_LB(object, type = "response")

  if(is.null(colnames(x))) colnames(x) <- paste0('V', 1:ncol(x))
  N1 <- ceiling(sum(x, na.rm = T))
  N0 <- length(as.matrix(x)) - N1

  TE <- lapply(seq(0, 1, length.out = ncuts), function(z){
    Pr <- ifelse(P>=z, 1, 0)
    c1 <- 1 - apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE)/apply(x == 1, 2, sum)
    c2 <- 1 - apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE)/apply(x == 0, 2, sum)
    TE <- 100/2 * (c1 + c2)
    lista <- list(TE)

  })
  TEp <-  data.frame(dplyr::bind_rows(TE), threshold = seq(0, 1, length.out = ncuts))

  thresholds <- TEp |>
                tidyr::pivot_longer(-threshold, names_to = "variable", values_to = "BACC") |>
                dplyr::group_by(variable) |>
                dplyr::mutate(merror = min(BACC)) |>
                dplyr::filter(BACC == merror) |>
                dplyr::mutate(row = dplyr::row_number()) |>
                dplyr::filter(row == 1) |> dplyr::ungroup() |>
                dplyr::select(variable, threshold, BACC)

  thresholds <- thresholds[match(colnames(x), thresholds$variable),]
  Pr <- matrix(NA, nrow(P), ncol(P))
  for(p in 1:ncol(P)){
    Pr[,p] <- ifelse( P[,p] >= thresholds$threshold[p], 1, 0)
  }
  c1 <- 1 - sum(apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE))/N1
  c2 <- 1 - sum(apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE))/N0
  BACC <- round(100/2 * (c1 + c2), 2)

  PCC = ifelse((x==1 & Pr==1) | (x==0 & Pr==0), 1, 0)
  PCC[is.na(PCC)] <- 0
  ones <-  apply(x, 2, sum, na.rm=TRUE)
  zeros <- length(as.matrix(x)) - ones

  confusion <- data.frame( Sensitivy = round(100*apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE)/ ones, 1),
                           Specificity = round(100*apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE)/ zeros, 1),
                           Global = round(100*colSums(PCC)/nrow(PCC), 1))

  out <- list(thresholds = thresholds, predictX = Pr, fitted = confusion, BER = BACC)
  class(out) <- c("BiplotML", "list")
  return(out)
}

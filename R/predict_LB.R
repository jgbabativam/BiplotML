#' @export
#' @title Predict Binary Responses from a Logistic Biplot
#' @description
#' Predicts the binary response matrix from a fitted logistic biplot and
#' computes the optimal classification threshold for each variable by minimising
#' the Balanced Error Rate (BER).
#'
#' @details
#' The optimal threshold for variable \eqn{j} is the value \eqn{\alpha_j \in [0,1]}
#' that minimises the Balanced Error Rate:
#' \deqn{BER_j = 1 - \frac{1}{2}
#'   \left(\frac{TP_j}{TP_j + FN_j} + \frac{TN_j}{TN_j + FP_j}\right),}
#' where \eqn{TP}, \eqn{TN}, \eqn{FP}, and \eqn{FN} denote true positives,
#' true negatives, false positives, and false negatives, respectively.
#'
#' @author Giovany Babativa <gbabativam@@gmail.com>
#'
#' @param object An object of class \code{BiplotML}, as returned by
#'   \code{\link{LogBip}}.
#' @param x The original binary matrix used to fit the model.
#' @param ncuts Number of equally spaced threshold candidates in \eqn{[0, 1]}.
#'   Default is \code{100}.
#'
#' @return A named list of class \code{BiplotML} with components:
#'   \describe{
#'     \item{\code{thresholds}}{Data frame with the optimal threshold and minimum
#'       BER for each variable.}
#'     \item{\code{predictX}}{Predicted binary matrix.}
#'     \item{\code{fitted}}{Confusion matrix (sensitivity, specificity, global
#'       accuracy) for each variable.}
#'     \item{\code{BER}}{Overall Balanced Error Rate (in percent).}
#'   }
#'
#' @examples
#' \donttest{
#' data("Methylation")
#' LB  <- LogBip(Methylation, plot = FALSE)
#' out <- pred_LB(LB, Methylation)
#' }

pred_LB <- function(object, x, ncuts = 100) {

  for (pkg in c("tidyr", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for this function. Please install it.",
           call. = FALSE)
    }
  }

  P <- fitted_LB(object, type = "response")

  if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))

  N1 <- ceiling(sum(x, na.rm = TRUE))
  N0 <- length(as.matrix(x)) - N1

  TE <- lapply(seq(0, 1, length.out = ncuts), function(z) {
    Pr <- ifelse(P >= z, 1, 0)
    c1 <- 1 - apply((Pr == 1) & (x == 1), 2, sum, na.rm = TRUE) /
              apply(x == 1, 2, sum)
    c2 <- 1 - apply((Pr == 0) & (x == 0), 2, sum, na.rm = TRUE) /
              apply(x == 0, 2, sum)
    list(100 / 2 * (c1 + c2))
  })

  TEp <- data.frame(dplyr::bind_rows(TE),
                    threshold = seq(0, 1, length.out = ncuts))

  thresholds <- TEp |>
    tidyr::pivot_longer(-threshold, names_to = "variable",
                        values_to = "BACC") |>
    dplyr::group_by(variable) |>
    dplyr::mutate(merror = min(BACC)) |>
    dplyr::filter(BACC == merror) |>
    dplyr::mutate(row = dplyr::row_number()) |>
    dplyr::filter(row == 1) |>
    dplyr::ungroup() |>
    dplyr::select(variable, threshold, BACC)

  thresholds <- thresholds[match(colnames(x), thresholds$variable), ]

  Pr <- matrix(NA, nrow(P), ncol(P))
  for (p in seq_len(ncol(P))) {
    Pr[, p] <- ifelse(P[, p] >= thresholds$threshold[p], 1, 0)
  }

  c1   <- 1 - sum(apply((Pr == 1) & (x == 1), 2, sum, na.rm = TRUE)) / N1
  c2   <- 1 - sum(apply((Pr == 0) & (x == 0), 2, sum, na.rm = TRUE)) / N0
  BACC <- round(100 / 2 * (c1 + c2), 2)

  PCC      <- ifelse((x == 1 & Pr == 1) | (x == 0 & Pr == 0), 1, 0)
  PCC[is.na(PCC)] <- 0
  ones     <- apply(x, 2, sum, na.rm = TRUE)
  zeros    <- length(as.matrix(x)) - ones

  confusion <- data.frame(
    Sensitivity = round(100 * apply((Pr == 1) & (x == 1), 2, sum,
                                    na.rm = TRUE) / ones,  1),
    Specificity = round(100 * apply((Pr == 0) & (x == 0), 2, sum,
                                    na.rm = TRUE) / zeros, 1),
    Global      = round(100 * colSums(PCC) / nrow(PCC), 1)
  )

  out <- list(thresholds = thresholds, predictX = Pr,
              fitted = confusion, BER = BACC)
  class(out) <- c("BiplotML", "list")
  return(out)
}

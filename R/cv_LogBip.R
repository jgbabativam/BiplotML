#' @export
#' @title Cross-Validation for Logistic Biplot
#' @description
#' Performs k-fold cross-validation for a logistic biplot model across a range
#' of dimensions, enabling selection of the optimal number of latent dimensions.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param data A binary matrix.
#' @param k Integer vector of dimensions to evaluate. Default is \code{0:5}.
#' @param K Number of folds. Default is \code{K = 7}.
#' @param method Fitting algorithm: \code{"MM"} (default), \code{"CG"},
#'   \code{"PDLB"}, or \code{"BFGS"}.
#' @param type Update formula for the CG method (see \code{\link{LogBip}}).
#' @param plot Logical; if \code{TRUE} (default), the cross-validation error
#'   curve is plotted.
#' @param maxit Maximum number of iterations. Defaults to \code{100} for
#'   gradient methods and \code{2000} for the MM algorithm.
#'
#' @return A data frame with columns \code{k}, \code{cv-error} (mean
#'   cross-validation error, in percent), and \code{train-error} (mean training
#'   error, in percent).
#'
#' @references
#' Bro, R., Kjeldahl, K., & Smilde, A. K. (2008). Cross-validation of component
#' models: a critical look at current methods. \emph{Analytical and Bioanalytical
#' Chemistry}, \emph{390}(5), 1241--1251.
#'
#' Wold, S. (1978). Cross-validatory estimation of the number of components in
#' factor and principal components models. \emph{Technometrics}, \emph{20}(4),
#' 397--405.
#'
#' @seealso \code{\link{LogBip}}, \code{\link{pred_LB}}, \code{\link{fitted_LB}},
#'   \code{\link{simBin}}
#'
#' @examples
#' \donttest{
#' set.seed(1234)
#' x <- simBin(n = 100, p = 50, k = 3, D = 0.5, C = 20)
#'
#' # Cross-validation using the MM algorithm
#' cv_MM <- cv_LogBip(data = x$X, k = 0:5, method = "MM", maxit = 1000)
#'
#' # Cross-validation using the Fletcher-Reeves CG algorithm
#' cv_CG <- cv_LogBip(data = x$X, k = 0:5, method = "CG", type = 1)
#'
#' # Cross-validation using the PDLB algorithm
#' cv_PB <- cv_LogBip(data = x$X, k = 0:5, method = "PDLB", maxit = 1000)
#' }

cv_LogBip <- function(data, k = 0:5, K = 7, method = "MM", type = NULL,
                      plot = TRUE, maxit = NULL) {

  x <- as.matrix(data)

  verify <- apply(x, 2, sd, na.rm = TRUE)
  if (any(verify == 0)) {
    stop("Some variables have zero variance; the procedure cannot be applied.")
  }

  min_k <- min(k)
  k_seq <- k

  if (any(is.na(x)) & method != "PDLB") {
    warning("The binary matrix contains missing values; method '", method,
            "' has been switched to 'PDLB'.")
    method <- "PDLB"
  }

  cvD <- matrix(NA, length(k_seq), 3)

  for (k in k_seq) {
    if (method == "CG" & k > 0) {
      bip   <- BiplotML::LogBip(x, k = k, method = method, type = type,
                                plot = FALSE)
      thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds

    } else if (method == "BFGS" & k > 0) {
      bip   <- BiplotML::LogBip(x, k = k, method = method, plot = FALSE)
      thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds

    } else if (method == "MM" & k > 0) {
      bip <- if (!is.null(maxit)) {
        BiplotML::LogBip(x, k = k, method = "MM", maxit = maxit, plot = FALSE)
      } else {
        BiplotML::LogBip(x, k = k, method = "MM", plot = FALSE)
      }
      thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds

    } else if (method == "PDLB" & k > 0) {
      bip <- if (!is.null(maxit)) {
        BiplotML::LogBip(x, k = k, method = "PDLB", maxit = maxit, plot = FALSE)
      } else {
        BiplotML::LogBip(x, k = k, method = "PDLB", plot = FALSE)
      }
      xt    <- bip$impute_x
      thres <- BiplotML::pred_LB(bip, bip$impute_x, ncuts = 50)$thresholds

    } else {
      theta <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(x, na.rm = TRUE)))
      P     <- plogis(theta)
      thres <- thresholds(x = x, P = P, ncuts = 50)$thres
    }

    if (k == 0 & any(is.na(x))) {
      m  <- apply(x, 2, mean, na.rm = TRUE)
      xt <- ifelse(is.na(x) & m > 0.5, 1, 0)
    }

    if (method != "PDLB")                        xt <- x
    if (k == 0 & any(!is.na(x)) & method == "PDLB") xt <- x

    folds    <- train_miss_pattern(xt, K = K)
    missVal  <- folds$missWold
    train    <- folds$Xtrain

    cv_errD  <- list()
    train_err <- list()

    for (i in seq_len(K)) {
      if (method == "CG" & k > 0) {
        missBip <- BiplotML::LogBip(train[[i]], k = k, method = method,
                                    type = type, plot = FALSE,
                                    cv_LogBip = TRUE)
        P <- fitted_LB(missBip, type = "response")

      } else if (method %in% c("BFGS", "MM", "PDLB") & k > 0) {
        missBip <- if (!is.null(maxit)) {
          BiplotML::LogBip(train[[i]], k = k, method = method, maxit = maxit,
                           plot = FALSE, cv_LogBip = TRUE)
        } else {
          BiplotML::LogBip(train[[i]], k = k, method = method, plot = FALSE,
                           cv_LogBip = TRUE)
        }
        P <- fitted_LB(missBip, type = "response")

      } else {
        theta <- rep(1, nrow(x)) %*%
          t(as.matrix(colMeans(train[[i]], na.rm = TRUE)))
        P <- plogis(theta)
      }

      Xhat <- matrix(NA, nrow(P), ncol(P))
      for (p in seq_len(ncol(P))) {
        Xhat[, p] <- ifelse(P[, p] >= thres$threshold[p], 1, 0)
      }
      Xhat[is.na(Xhat)] <- 0

      Xhat_pred <- Xhat[missVal[[i]]]
      xReal     <- xt[missVal[[i]]]

      n1  <- sum(xReal);           n1t <- sum(xt)
      n0  <- length(xReal) - n1;   n0t <- length(xt) - n1t

      err0 <- sum(ifelse(xReal == 0 & Xhat_pred == 1, 1, 0)) / n0
      err1 <- sum(ifelse(xReal == 1 & Xhat_pred == 0, 1, 0)) / n1

      err0t <- sum(ifelse(xt == 0 & Xhat == 1, 1, 0)) / n0t
      err1t <- sum(ifelse(xt == 1 & Xhat == 0, 1, 0)) / n1t

      cv_errD[[i]]   <- 100 / 2 * (err0 + err1)
      train_err[[i]] <- 100 / 2 * (err0t + err1t)
    }

    row <- k - min_k + 1
    cvD[row, 1] <- k
    cvD[row, 2] <- round(mean(sapply(cv_errD,   mean, na.rm = TRUE),
                              na.rm = TRUE), 2)
    cvD[row, 3] <- round(mean(sapply(train_err, mean, na.rm = TRUE),
                              na.rm = TRUE), 2)
  }

  out <- as.data.frame(cvD)
  colnames(out) <- c("k", "cv-error", "train-error")

  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE) |
        !requireNamespace("tidyr",   quietly = TRUE)) {
      message("Packages 'ggplot2' and 'tidyr' are required for plotting.")
    } else {
      cvp <- out |>
        tidyr::pivot_longer(-k, names_to = "type.error", values_to = "error") |>
        ggplot2::ggplot(ggplot2::aes(x = k, y = error,
                                     linetype = type.error,
                                     shape    = type.error,
                                     color    = type.error)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 3.0) +
        ggplot2::geom_vline(
          xintercept = out[which.min(out$`cv-error`), "k"],
          linetype = 2) +
        ggplot2::scale_color_manual(values = c("red", "blue")) +
        ggplot2::labs(x = "Dimensions (k)", y = "Classification error (%)",
                      caption = paste("Cross-validation using the",
                                      method, "method")) +
        ggplot2::scale_x_continuous(breaks = seq(0, 20, 1)) +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position  = "top",
                       legend.title     = ggplot2::element_blank(),
                       legend.text      = ggplot2::element_text(size = 11))
      print(cvp)
    }
  }

  return(out)
}

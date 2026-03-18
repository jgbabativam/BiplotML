#' @title Fitted Values for a Logistic Biplot
#' @description
#' Computes the fitted (predicted) matrix for a logistic biplot model on either
#' the logit (log-odds) scale or the probability scale.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param object An object of class \code{BiplotML}, as returned by
#'   \code{\link{LogBip}}.
#' @param type Scale of the fitted values: \code{"link"} for the logit scale
#'   (log-odds) or \code{"response"} for the probability scale. Partial matching
#'   is supported.
#'
#' @return A numeric matrix of fitted values with the same dimensions as the
#'   original binary matrix.
#'
#' @examples
#' \donttest{
#' data("Methylation")
#' LB    <- LogBip(Methylation, plot = FALSE)
#' Theta <- fitted_LB(LB, type = "link")      # log-odds scale
#' Pi    <- fitted_LB(LB, type = "response")  # probability scale
#' }
#' @export fitted_LB

fitted_LB <- function(object, type = c("link", "response")) {
  type  <- match.arg(type[1], c("link", "response"))
  n     <- nrow(object$Ahat)
  theta <- as.matrix(cbind(rep(1, n), object$Ahat)) %*% t(object$Bhat)
  P     <- plogis(theta)

  if (type == "link") {
    return(theta)
  } else {
    return(P)
  }
}

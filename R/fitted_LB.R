#' @title Fitted values using Logistic Biplot
#'
#' @description
#' Fit a lower dimentional representation of the binary matrix using logistic biplot
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param object BiplotML object
#' @param type the type of fitting required. \code{type = "link"} gives output on the logit scale and
#'  \code{type = "response"} gives output on the probability scale
#' @examples
#' \dontrun{
#' data("Methylation")
#' LB <- LogBip(Methylation, plot = FALSE)
#' Theta <- fitted_LB(LB, type = "link")
#' Pi <- fitted_LB(LB, type = "response")
#' }
#' @export fitted_LB

fitted_LB <- function(object, type = c("link", "response")){
  type <- match.arg(type[1], c("link","response"))
  n = nrow(object$Ahat)

  theta = as.matrix(cbind(rep(1,n),object$Ahat))%*%t(object$Bhat)
  P = plogis(theta)

  if (type == "link") {
    return(theta)
  } else if (type == "response") {
    return(P)
  }
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted_LB.R
\name{fitted_LB}
\alias{fitted_LB}
\title{Fitted values using Logistic Biplot}
\usage{
fitted_LB(object, type = c("link", "response"))
}
\arguments{
\item{object}{BiplotML object}

\item{type}{the type of fitting required. \code{type = "link"} gives output on the logit scale and
\code{type = "response"} gives output on the probability scale}
}
\description{
Fit a lower-dimensional representation of the binary matrix using logistic biplot
}
\examples{
\dontrun{
data("Methylation")
LB <- LogBip(Methylation, plot = FALSE)
Theta <- fitted_LB(LB, type = "link")
Pi <- fitted_LB(LB, type = "response")
}
}
\author{
Giovany Babativa <gbabativam@gmail.com>
}

#' DNA Methylation Binary Data
#'
#' A binary matrix of DNA methylation measurements for a sample of individuals.
#' Each row represents an individual and each column a CpG site; a value of 1
#' indicates methylation and 0 indicates no methylation.
#'
#' @format A binary matrix with 50 rows (individuals) and 13 columns (CpG sites).
#'
#' @source Publicly available methylation data used for illustrative purposes.
#'
#' @examples
#' data("Methylation")
#' dim(Methylation)
"Methylation"

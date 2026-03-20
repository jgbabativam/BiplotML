#' @importFrom stats IQR
#'
#' @export
#' @title Plot a Binary Logistic Biplot
#' @description
#' Produces a ggplot2-based logistic biplot from a \code{BiplotML} object,
#' optionally drawing confidence ellipses for the row markers.
#'
#' @details
#' Variable vectors are drawn as arrows from the point where the predicted
#' probability equals 0.5 to the point where it equals \code{endsegm}.
#' Short arrows indicate a rapid increase in the probability of the
#' corresponding characteristic. The orthogonal projection of a row marker onto
#' a variable's arrow approximates the probability that the characteristic is
#' present for that individual.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param x An object of class \code{BiplotML}, as returned by
#'   \code{\link{LogBip}} or \code{\link{bootBLB}}.
#' @param dim Integer vector of length 2 specifying which dimensions to plot.
#'   Default is \code{c(1, 2)}.
#' @param col.ind Color for the row markers. If \code{NULL} (default), a single
#'   neutral color is used.
#' @param col.var Color for the variable arrows. Default is \code{"#0E185F"}.
#' @param label.ind Logical; if \code{TRUE}, row points are labelled. Default
#'   is \code{FALSE}.
#' @param draw Which graph to draw: \code{"biplot"} (default) for both rows and
#'   columns, \code{"ind"} for individuals only, or \code{"var"} for variables
#'   only.
#' @param titles Main title for the plot. Default is \code{NULL} (no title).
#' @param subtitle Subtitle for the plot, typically the algorithm name. Default
#'   is \code{NULL}.
#' @param ellipses Logical; if \code{TRUE}, confidence ellipses are drawn around
#'   the row markers (requires a bootstrap fit from \code{\link{bootBLB}}).
#'   Default is \code{FALSE}.
#' @param endsegm End point of the variable arrow on the probability scale.
#'   The arrow starts at 0.5 and ends at this value. Default is \code{0.75}.
#' @param repel Logical; if \code{TRUE}, overlapping row labels are repelled
#'   using \pkg{ggrepel}. Default is \code{FALSE}.
#' @param xylim Numeric vector of length 2 specifying the common range of both
#'   axes, e.g. \code{c(-10, 10)}. If \code{NULL} (default), the range is
#'   determined automatically from all plotted elements including ellipses.
#'
#' @return A \code{ggplot2} object.
#'
#' @references
#' Meulman, J. J., & Heiser, W. J. (1983). \emph{The Display of Bootstrap
#' Solutions in Multidimensional Scaling} (Technical memorandum). Bell
#' Laboratories.
#'
#' Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
#' M. Greenacre & J. Blasius (Eds.), \emph{Multiple Correspondence Analysis and
#' Related Methods} (pp. 503--521). Chapman & Hall.
#'
#' @seealso \code{\link{bootBLB}}
#'
#' @examples
#' \donttest{
#' data("Methylation")
#' set.seed(123456)
#' res <- LogBip(x = Methylation, method = "MM", maxit = 1000, plot = FALSE)
#' plotBLB(x = res, titles = "Methylation Logistic Biplot", endsegm = 0.95)
#' }



plotBLB <- function(x, dim = c(1, 2), col.ind = NULL, col.var = "#0E185F",
                    label.ind = FALSE,
                    draw     = c("biplot", "ind", "var"),
                    titles   = NULL, subtitle = NULL,
                    ellipses = FALSE,
                    endsegm  = 0.75, repel = FALSE, xylim = NULL) {

  for (pkg in c("ggplot2", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for this function. Please install it.",
           call. = FALSE)
    }
  }

  EspA <- x$Ahat
  EspB <- x$Bhat
  k    <- ncol(EspA)

  if (k < max(dim)) {
    stop("The logistic biplot has only ", k,
         " dimension(s). Reduce the requested plot dimensions.")
  }

  grap <- match.arg(draw[1], c("ind", "var", "biplot"))

  d1  <- paste0("Dim", dim[1])
  d2  <- paste0("Dim", dim[2])
  dd1 <- paste0("bb",  dim[1])
  dd2 <- paste0("bb",  dim[2])

  denom <- rowSums(EspB[, c(dd1, dd2), drop = FALSE] *
                     EspB[, c(dd1, dd2), drop = FALSE])

  lp <- log(endsegm / (1 - endsegm))

  EspB[["x.50"]]  <- (-EspB[["bb0"]] * EspB[[dd1]]) / denom
  EspB[["y.50"]]  <- (-EspB[["bb0"]] * EspB[[dd2]]) / denom
  EspB[["x.end"]] <- (lp - EspB[["bb0"]]) * EspB[[dd1]] / denom
  EspB[["y.end"]] <- (lp - EspB[["bb0"]]) * EspB[[dd2]] / denom

  EspA[["label"]] <- rownames(EspA)
  EspB[["label"]] <- rownames(EspB)

  # Compute axis limits including ellipse coordinates if present,
  # then clip extreme ellipse outliers to 3 * IQR beyond the main data range
  # to prevent a few runaway ellipses from collapsing the rest of the plot.
  if (!is.null(xylim)) {
    lims <- xylim
  } else {
    core_vals <- c(EspA[[d1]], EspA[[d2]],
                   EspB[["x.50"]],  EspB[["y.50"]],
                   EspB[["x.end"]], EspB[["y.end"]])
    core_range <- range(core_vals, na.rm = TRUE)
    iqr_val    <- IQR(core_vals, na.rm = TRUE)
    # Allow ellipses to extend up to 3 IQR beyond the core range
    expand     <- max(3 * iqr_val, diff(core_range) * 0.15)
    lims       <- core_range + c(-expand, expand)

    # If ellipses are being drawn, include their coords but cap at the
    # expanded limits so outlier ellipses do not blow up the scale
    if (ellipses && !is.null(x$Ellip)) {
      ellip_vals <- c(x$Ellip[[d1]], x$Ellip[[d2]])
      ellip_range <- range(ellip_vals, na.rm = TRUE)
      lims[1] <- max(lims[1], ellip_range[1])
      lims[2] <- min(lims[2], ellip_range[2])
      # Safety: ensure lims[1] < lims[2]
      if (lims[1] >= lims[2]) lims <- core_range + c(-expand, expand)
    }
    margin <- diff(lims) * 0.04
    lims   <- lims + c(-margin, margin)
  }

  col_ind <- if (is.null(col.ind)) "#444444" else col.ind

  # Build title: use subtitle if provided, else fall back to x$method
  auto_subtitle <- if (!is.null(subtitle)) subtitle else x$method
  plot_title    <- if (!is.null(auto_subtitle) && !is.null(titles)) {
    paste0(titles, "\n", auto_subtitle)
  } else if (!is.null(titles)) {
    titles
  } else {
    auto_subtitle
  }

  # Base plot
  p_base <- ggplot2::ggplot() +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(xlim = lims, ylim = lims) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey60") +
    ggplot2::labs(title    = plot_title,
                  x = paste0("Dimension ", dim[1]),
                  y = paste0("Dimension ", dim[2])) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      hjust = 0.5, size = 11, face = "bold"))

  # Confidence ellipses (drawn FIRST so they appear behind the points)
  if (ellipses && !is.null(x$Ellip)) {
    Ellip_data <- x$Ellip
    p_base <- p_base +
      ggplot2::geom_path(
        data = Ellip_data,
        ggplot2::aes(x = .data[[d1]], y = .data[[d2]],
                     group = .data[["ind"]]),
        color = col_ind, linewidth = 0.3, alpha = 0.5)
  }

  # Row markers
  if (grap %in% c("ind", "biplot")) {
    p_base <- p_base +
      ggplot2::geom_point(data = EspA,
                          ggplot2::aes(x = .data[[d1]],
                                       y = .data[[d2]]),
                          color = col_ind, size = 1.5)
    if (label.ind) {
      if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
        p_base <- p_base +
          ggrepel::geom_text_repel(
            data = EspA,
            ggplot2::aes(x = .data[[d1]], y = .data[[d2]],
                         label = .data[["label"]]),
            size = 3, color = col_ind)
      } else {
        p_base <- p_base +
          ggplot2::geom_text(
            data = EspA,
            ggplot2::aes(x = .data[[d1]], y = .data[[d2]],
                         label = .data[["label"]]),
            size = 3, color = col_ind, vjust = -0.5)
      }
    }
  }

  # Variable arrows
  if (grap %in% c("var", "biplot")) {
    p_base <- p_base +
      ggplot2::geom_segment(
        data = EspB,
        ggplot2::aes(x    = .data[["x.50"]],  y    = .data[["y.50"]],
                     xend = .data[["x.end"]], yend = .data[["y.end"]]),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
        color = col.var) +
      ggplot2::geom_text(
        data = EspB,
        ggplot2::aes(x     = .data[["x.end"]],
                     y     = .data[["y.end"]],
                     label = .data[["label"]]),
        size = 3, color = col.var, vjust = -0.4)
  }

  return(p_base)
}

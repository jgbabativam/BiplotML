#' @importFrom stats IQR
#'
#' @export
#' @title Plot a Binary Logistic Biplot
#' @description
#' Produces a \pkg{ggplot2}-based logistic biplot from a \code{BiplotML}
#' object fitted with \code{\link{LogBip}}. Supports coloring and
#' shaping of row markers by a categorical variable, filled arrowheads,
#' dashed reference lines that span the full plot area, and flexible axis-limit control via
#' \code{xylim}, \code{xlim}, and \code{ylim}.
#'
#' @details
#' Variable vectors are drawn as arrows from the point where the predicted
#' probability equals 0.5 to the point where it equals \code{endsegm}.
#' Short arrows indicate a rapid increase in the probability of the
#' corresponding characteristic. The orthogonal projection of a row marker
#' onto a variable's arrow approximates the probability that the
#' characteristic is present for that individual.
#'
#' The three arguments that control axis limits are evaluated in the
#' following order of priority:
#' \enumerate{
#'   \item \code{xlim} and \code{ylim} (independent limits for each axis).
#'   \item \code{xylim} (symmetric limits applied to both axes).
#'   \item Automatic limits derived from all plotted elements.
#' }
#'
#' The \code{escala} argument multiplies the row marker coordinates before
#' plotting so that they are visually comparable to the variable arrows,
#' which are expressed in the original parameter units. It only affects the
#' display, not the stored coordinates.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param x An object of class \code{BiplotML}, as returned by
#'   \code{\link{LogBip}}.
#' @param dim Integer vector of length 2 specifying which dimensions to plot.
#'   Default is \code{c(1, 2)}.
#' @param col.ind Optional vector of the same length as the number of rows in
#'   the original data, used to color \emph{and} shape the row markers by a
#'   categorical variable (e.g., \code{col.ind = df$group}). Levels are
#'   mapped to the \code{"Set1"} palette and to filled geometric shapes.
#'   If \code{NULL} (default), all row markers are drawn as gold triangles
#'   (\code{shape = 17}, color \code{"#E7B800"}) when no \code{col.ind} is
#'   provided.
#' @param col.var Color for the variable arrows. Default is \code{"#0E185F"}
#'   (dark navy).
#' @param label.ind Logical; if \code{TRUE}, row markers are labelled.
#'   Default is \code{FALSE}.
#' @param draw Which graph to draw. One of \code{"biplot"} (default, both row
#'   and column markers), \code{"ind"} (row markers only), or \code{"var"}
#'   (variable arrows only). Partial matching is supported.
#' @param titles Main title for the plot. If \code{NULL} (default), a
#'   generic title is used depending on \code{draw}.
#' @param ellipses Logical; if \code{TRUE}, bootstrap confidence ellipses are
#'   drawn around the row markers. Requires a bootstrap fit
#' @param endsegm End point of the variable arrow on the probability scale.
#'   The arrow starts at \eqn{p = 0.5} and ends at this value.
#'   Default is \code{0.75}.
#' @param repel Logical; if \code{TRUE}, overlapping variable labels are
#'   repelled using \pkg{ggrepel}. Default is \code{FALSE}.
#' @param xylim Numeric vector of length 2 specifying a symmetric range
#'   applied to both axes, e.g., \code{c(-80, 80)}. Overrides automatic
#'   limits. Takes precedence over automatic limits but is overridden by
#'   \code{xlim}/\code{ylim} if those are also supplied.
#'   Default is \code{NULL}.
#' @param xlim Numeric vector of length 2 specifying the range of the
#'   x-axis independently, e.g., \code{c(-100, 60)}. Takes precedence over
#'   \code{xylim}. Default is \code{NULL}.
#' @param ylim Numeric vector of length 2 specifying the range of the
#'   y-axis independently, e.g., \code{c(-80, 80)}. Takes precedence over
#'   \code{xylim}. Default is \code{NULL}.
#' @param escala Positive numeric scalar. Multiplicative factor applied to
#'   the row marker coordinates (\code{x$Ahat}) before plotting, so that
#'   they are on a comparable visual scale to the variable arrows.
#'   If \code{NULL} (default), the value is chosen automatically so that
#'   the range of the scaled row markers matches the range of the variable
#'   arrows, producing a visually balanced biplot. Pass an explicit numeric
#'   value to override the automatic calculation (e.g., \code{escala = 65}).
#'
#' @return A \code{ggplot2} object that can be further customised with
#'   standard \pkg{ggplot2} functions (e.g., \code{theme()}, \code{labs()}).
#'
#' @references
#' Meulman, J. J., & Heiser, W. J. (1983). \emph{The Display of Bootstrap
#' Solutions in Multidimensional Scaling} (Technical memorandum). Bell
#' Laboratories.
#'
#' Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
#' M. Greenacre & J. Blasius (Eds.), \emph{Multiple Correspondence Analysis
#' and Related Methods} (pp. 503--521). Chapman & Hall.
#'
#' @seealso \code{\link{LogBip}}
#'
#' @examples
#' \donttest{
#' data("Methylation")
#' set.seed(123456)
#' res <- LogBip(x = Methylation, method = "MM", maxit = 1000, plot = FALSE)
#' }

plotBLB <- function(x,
                    dim      = c(1, 2),
                    col.ind  = NULL,
                    col.var  = "#0E185F",
                    label.ind = FALSE,
                    draw     = c("biplot", "ind", "var"),
                    titles   = NULL,
                    ellipses = FALSE,
                    endsegm  = 0.75,
                    repel    = FALSE,
                    xylim    = NULL,
                    xlim     = NULL,
                    ylim     = NULL,
                    escala   = NULL) {

  for (pkg in c("ggplot2", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Please install it.",
           call. = FALSE)
    }
  }

  # ── Extract fitting objects ────────────────────────────────────────────────
  EspA <- x$Ahat
  EspB <- x$Bhat
  k    <- ncol(EspA)

  if (k < max(dim)) {
    stop("The logistic biplot has only ", k,
         " dimension(s). Reduce the requested plot dimensions.")
  }

  grap <- match.arg(draw[1], c("ind", "var", "biplot"))

  # ── Dimension column names ─────────────────────────────────────────────────
  colnames(EspA) <- paste0("Dim", seq_len(k))
  d1  <- paste0("Dim", dim[1])
  d2  <- paste0("Dim", dim[2])
  dd1 <- paste0("bb",  dim[1])
  dd2 <- paste0("bb",  dim[2])

  # ── Column marker (variable) coordinates ──────────────────────────────────
  # Arrow starts at p = 0.5 and ends at p = endsegm
  denom <- rowSums(EspB[, c(dd1, dd2), drop = FALSE] *
                     EspB[, c(dd1, dd2), drop = FALSE])
  lp <- log(endsegm / (1 - endsegm))

  EspB[["x.50"]]  <- (-EspB[["bb0"]] * EspB[[dd1]]) / denom
  EspB[["y.50"]]  <- (-EspB[["bb0"]] * EspB[[dd2]]) / denom
  EspB[["x.end"]] <- (lp - EspB[["bb0"]]) * EspB[[dd1]] / denom
  EspB[["y.end"]] <- (lp - EspB[["bb0"]]) * EspB[[dd2]] / denom

  EspB[["label"]] <- rownames(EspB)
  EspA[["label"]] <- rownames(EspA)

  # ── Row marker scaled coordinates ─────────────────────────────────────────
  if (is.null(escala)) {
    range_ind <- max(
      diff(range(c(EspA[[d1]], EspA[[d2]]), na.rm = TRUE)),
      1e-6
    )
    range_var <- diff(range(
      c(EspB[["x.50"]], EspB[["x.end"]],
        EspB[["y.50"]], EspB[["y.end"]]),
      na.rm = TRUE
    ))
    escala <- range_var / range_ind
  }

  EspA[["x.sc"]] <- EspA[[d1]] * escala
  EspA[["y.sc"]] <- EspA[[d2]] * escala

  # ── Axis limits ────────────────────────────────────────────────────────────
  # Priority: xlim/ylim > xylim > automatic
  if (!is.null(xlim) && !is.null(ylim)) {
    x.lim <- xlim
    y.lim <- ylim
  } else if (!is.null(xylim)) {
    x.lim <- xylim
    y.lim <- xylim
  } else {
    # Automatic: derive from all plotted elements
    if (ellipses && !is.null(x$Ellip)) {
      all_x <- c(EspA[["x.sc"]], EspB[["x.50"]], EspB[["x.end"]],
                 x$Ellip[[d1]])
      all_y <- c(EspA[["y.sc"]], EspB[["y.50"]], EspB[["y.end"]],
                 x$Ellip[[d2]])
    } else {
      all_x <- c(EspA[["x.sc"]], EspB[["x.50"]], EspB[["x.end"]])
      all_y <- c(EspA[["y.sc"]], EspB[["y.50"]], EspB[["y.end"]])
    }
    margin <- 0.05
    x.lim  <- range(all_x, na.rm = TRUE) * c(1 + margin, 1 + margin)
    y.lim  <- range(all_y, na.rm = TRUE) * c(1 + margin, 1 + margin)
    # Symmetric auto limits: use the wider range for both axes so
    # coord_fixed() does not distort the biplot geometry
    sym    <- max(abs(c(x.lim, y.lim)))
    x.lim  <- c(-sym, sym)
    y.lim  <- c(-sym, sym)
  }

  # ── Subtitle from method ───────────────────────────────────────────────────
  if (x$method == "CG") {
    subt <- "Estimation with Conjugate Gradient"
  } else {
    subt <- paste0("Estimation with ", x$method, " algorithm")
  }

  # ── Default titles per draw mode ───────────────────────────────────────────
  if (is.null(titles)) {
    titulo <- switch(grap,
                     ind    = "Individuals plot",
                     var    = "Variables plot",
                     biplot = "Logistic Biplot")
  } else {
    titulo <- titles
  }

  make_arrow <- function(angle = 20, len = 0.3) {
    grid::arrow(angle  = angle,
                type   = "closed",
                ends   = "last",
                length = grid::unit(len, "cm"))
  }

  ref_lines <- list(
    ggplot2::geom_vline(xintercept = 0, color = "black",
                        linetype = "dashed", linewidth = 0.4),
    ggplot2::geom_hline(yintercept = 0, color = "black",
                        linetype = "dashed", linewidth = 0.4)
  )

  g <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ref_lines +
    ggplot2::xlim(x.lim[1], x.lim[2]) +
    ggplot2::ylim(y.lim[1], y.lim[2]) +
    ggplot2::xlab(paste("Dimension", dim[1])) +
    ggplot2::ylab(paste("Dimension", dim[2])) +
    ggplot2::labs(title = titulo, subtitle = subt) +
    ggplot2::theme(
      axis.text     = ggplot2::element_text(face = "bold"),
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  # ══════════════════════════════════════════════════════════════════════════
  # IND plot
  # ══════════════════════════════════════════════════════════════════════════
  if (grap == "ind") {

    if (is.null(col.ind)) {
      g <- g +
        ggplot2::geom_point(
          data  = EspA,
          ggplot2::aes(x = .data[["x.sc"]], y = .data[["y.sc"]]),
          color = "#E7B800", size = 2, shape = 17
        ) +
        ggplot2::theme(legend.position = "none")
    } else {
      g <- g +
        ggplot2::geom_point(
          data  = EspA,
          ggplot2::aes(x      = .data[["x.sc"]],
                       y      = .data[["y.sc"]],
                       colour = col.ind,
                       shape  = col.ind),
          size = 2
        ) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::scale_shape_manual(
          values = c(16, 17, 15, 18, 19, 3, 4, 8)[
            seq_along(levels(factor(col.ind)))]
        )
    }

    if (label.ind) {
      if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
        g <- g + ggrepel::geom_text_repel(
          data  = EspA,
          ggplot2::aes(x     = .data[["x.sc"]],
                       y     = .data[["y.sc"]],
                       label = .data[["label"]]),
          size          = 3.5,
          segment.color = "grey50"
        )
      } else {
        g <- g + ggplot2::geom_text(
          data  = EspA,
          ggplot2::aes(x     = .data[["x.sc"]],
                       y     = .data[["y.sc"]],
                       label = .data[["label"]]),
          size = 3.5
        )
      }
    }
  }

  # ══════════════════════════════════════════════════════════════════════════
  # VAR plot
  # ══════════════════════════════════════════════════════════════════════════
  if (grap == "var") {

    g <- g +
      ggplot2::geom_segment(
        data     = EspB,
        ggplot2::aes(x    = .data[["x.50"]],
                     y    = .data[["y.50"]],
                     xend = .data[["x.end"]],
                     yend = .data[["y.end"]]),
        arrow       = make_arrow(angle = 25),
        colour      = col.var,
        lineend     = "butt",
        linejoin    = "mitre",
        show.legend = FALSE
      )

    if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
      g <- g + ggrepel::geom_text_repel(
        data  = EspB,
        ggplot2::aes(x     = .data[["x.end"]],
                     y     = .data[["y.end"]],
                     label = .data[["label"]],
                     angle = (180 / pi) * atan(.data[["y.end"]] /
                                                 .data[["x.end"]]),
                     hjust = (1 - 2 * sign(.data[["x.end"]])) / 2),
        size          = 3,
        segment.color = "grey50",
        colour        = "black",
        vjust         = 0
      )
    } else {
      g <- g + ggplot2::geom_text(
        data  = EspB,
        ggplot2::aes(x     = .data[["x.end"]],
                     y     = .data[["y.end"]],
                     label = .data[["label"]],
                     angle = (180 / pi) * atan(.data[["y.end"]] /
                                                 .data[["x.end"]]),
                     hjust = (1 - 2 * sign(.data[["x.end"]])) / 2),
        size   = 3,
        colour = "black",
        vjust  = 0
      )
    }
  }

  # ══════════════════════════════════════════════════════════════════════════
  # BIPLOT
  # ══════════════════════════════════════════════════════════════════════════
  if (grap == "biplot") {

    # ── Row markers ───────────────────────────────────────────────────────
    if (is.null(col.ind)) {
      g <- g +
        ggplot2::geom_point(
          data  = EspA,
          ggplot2::aes(x = .data[["x.sc"]], y = .data[["y.sc"]]),
          color = "#E7B800", size = 2, shape = 17
        ) +
        ggplot2::theme(legend.position = "none")
    } else {
      g <- g +
        ggplot2::geom_point(
          data  = EspA,
          ggplot2::aes(x      = .data[["x.sc"]],
                       y      = .data[["y.sc"]],
                       colour = col.ind,
                       shape  = col.ind),
          size = 2
        ) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::scale_shape_manual(
          values = c(16, 17, 15, 18, 19, 3, 4, 8)[
            seq_along(levels(factor(col.ind)))]
        )
    }

    # ── Variable arrows ───────────────────────────────────────────────────
    g <- g +
      ggplot2::geom_segment(
        data     = EspB,
        ggplot2::aes(x    = .data[["x.50"]],
                     y    = .data[["y.50"]],
                     xend = .data[["x.end"]],
                     yend = .data[["y.end"]]),
        arrow       = make_arrow(angle = 20),
        colour      = col.var,
        lineend     = "butt",
        linejoin    = "mitre",
        show.legend = FALSE
      )

    # ── Variable labels ───────────────────────────────────────────────────
    if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
      g <- g + ggrepel::geom_text_repel(
        data  = EspB,
        ggplot2::aes(x     = .data[["x.end"]],
                     y     = .data[["y.end"]],
                     label = .data[["label"]],
                     angle = (180 / pi) * atan(.data[["y.end"]] /
                                                 .data[["x.end"]]),
                     hjust = (1 - 2 * sign(.data[["x.end"]])) / 2),
        size          = 3,
        segment.color = "grey50",
        vjust         = 0
      )
    } else {
      g <- g + ggplot2::geom_text(
        data  = EspB,
        ggplot2::aes(x     = .data[["x.end"]],
                     y     = .data[["y.end"]],
                     label = .data[["label"]],
                     angle = (180 / pi) * atan(.data[["y.end"]] /
                                                 .data[["x.end"]]),
                     hjust = (1 - 1.0 * sign(.data[["x.end"]])) / 2),
        size  = 3,
        vjust = 0
      )
    }

    # ── Row labels ────────────────────────────────────────────────────────
    if (label.ind) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        g <- g + ggrepel::geom_text_repel(
          data  = EspA,
          ggplot2::aes(x     = .data[["x.sc"]],
                       y     = .data[["y.sc"]],
                       label = .data[["label"]]),
          size          = 3.5,
          segment.color = "grey50"
        )
      } else {
        g <- g + ggplot2::geom_text(
          data  = EspA,
          ggplot2::aes(x     = .data[["x.sc"]],
                       y     = .data[["y.sc"]],
                       label = .data[["label"]]),
          size  = 3.5,
          vjust = -0.5
        )
      }
    }
  }

  # ══════════════════════════════════════════════════════════════════════════
  # Bootstrap ellipses
  # ══════════════════════════════════════════════════════════════════════════
  if (ellipses && !is.null(x$Ellip)) {
    g <- g +
      ggplot2::geom_point(
        data  = x$Ellip,
        ggplot2::aes(x = .data[[d1]], y = .data[[d2]],
                     group = .data[["ind"]]),
        size   = 0.001,
        colour = "lightgray",
        shape  = 20
      )
  }

  return(g + ggplot2::coord_fixed())
}

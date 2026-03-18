#' @importFrom optimx optimr
#' @importFrom shapes procOPA
#' @importFrom stats runif
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @export
#' @title Bootstrap Binary Logistic Biplot
#' @description
#' Fits a binary logistic biplot and applies a bootstrap procedure to obtain
#' confidence ellipses for the row markers.
#'
#' @details
#' When \code{sup = TRUE}, individuals not selected in a bootstrap resample are
#' treated as supplementary rows: their coordinates are estimated by fixing the
#' column markers \strong{B} from the resample fit and optimising only over the
#' row coordinates. When \code{sup = FALSE}, only the resampled individuals are
#' used and Procrustes rotation aligns bootstrap solutions to the full-data
#' solution.
#'
#' @author Giovany Babativa <jgbabativam@@unal.edu.co>
#'
#' @param x A binary matrix.
#' @param k Number of dimensions. Default is \code{k = 2}.
#' @param L Ridge penalization parameter. Default is \code{L = 0}.
#' @param method Fitting algorithm: \code{"CG"} (default), \code{"BFGS"}, or
#'   others accepted by \code{\link{LogBip}}.
#' @param type Update formula for the CG method (see \code{\link{LogBip}}).
#'   Default is \code{type = 1} (Fletcher-Reeves).
#' @param plot Logical; if \code{TRUE} (default), the bootstrap biplot is
#'   plotted after fitting.
#' @param sup Logical; if \code{TRUE} (default), out-of-sample individuals are
#'   projected as supplementary rows. See Details.
#' @param ellipses Logical; if \code{TRUE}, confidence ellipses are drawn around
#'   row markers. Default is \code{FALSE}.
#' @param maxit Maximum number of iterations.
#' @param resamples Number of bootstrap resamples. Default is \code{100}.
#' @param conf Confidence level for the ellipses. Default is \code{0.90}.
#' @param col.ind Color for the row markers. Passed to \code{\link{plotBLB}}.
#'
#' @return An object of class \code{BiplotML} (a named list) containing:
#'   \describe{
#'     \item{\code{Ahat}}{Mean row-marker coordinates across bootstrap resamples.}
#'     \item{\code{Bhat}}{Mean column-marker coordinates across bootstrap resamples.}
#'     \item{\code{BootA}}{Data frame of row-marker coordinates for all resamples.}
#'     \item{\code{BootB}}{Data frame of column-marker coordinates for all resamples.}
#'     \item{\code{pred}}{Predicted binary matrix.}
#'     \item{\code{fit}}{Confusion matrix (sensitivity, specificity, global accuracy).}
#'     \item{\code{rows}}{Row names.}
#'     \item{\code{cols}}{Column names.}
#'     \item{\code{method}}{Fitting method used.}
#'     \item{\code{Ellip}}{Data frame of ellipse coordinates (only when
#'       \code{ellipses = TRUE}).}
#'   }
#'
#' @references
#' Milan, L., & Whittaker, J. (1995). Application of the parametric bootstrap
#' to models that incorporate a singular value decomposition.
#' \emph{Applied Statistics}, \emph{44}, 31--49.
#'
#' Vicente-Villardon, J. L., & Galindo, M. P. (2006). Logistic biplots. In
#' M. Greenacre & J. Blasius (Eds.), \emph{Multiple Correspondence Analysis and
#' Related Methods} (pp. 503--521). Chapman & Hall.
#'
#' @seealso \code{\link{plotBLB}}, \code{\link{performanceBLB}}
#'
#' @examples
#' \donttest{
#' data("Methylation")
#' set.seed(02052020)
#' out_sup  <- bootBLB(x = Methylation, ellipses = FALSE)
#' out_boot <- bootBLB(x = Methylation, sup = FALSE, ellipses = TRUE)
#' }

bootBLB <- function(x, k = 2, L = 0, method = "CG", type = 1, plot = TRUE,
                    sup = TRUE, ellipses = FALSE, maxit = NULL,
                    resamples = 100, conf = 0.9, col.ind = NULL) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for this function. Please install it.",
         call. = FALSE)
  }

  n <- nrow(x); p <- ncol(x)
  aik <- n * k; bjk <- p * (k + 1)
  dTheta <- aik + bjk; s <- k + 1

  # Helper: run optimr suppressing the verbose CG control messages
  run_optimr <- function(par, xt, use_type = FALSE) {
    if (use_type) {
      suppressWarnings(
        optimr(par = par, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
               xt = xt, k = k, lambda = L, method = method,
               control = list(type = type))
      )
    } else {
      suppressWarnings(
        optimr(par = par, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
               xt = xt, k = k, lambda = L, method = method)
      )
    }
  }

  use_type <- (method == "CG")

  # --- Initial fit on full data ---
  res <- run_optimr(runif(dTheta), x, use_type)

  par <- res$par
  A   <- data.frame(matrix(res$par[1:aik], n, k))
  colnames(A)  <- paste0("Dim", seq_len(k))
  rownames(A)  <- rownames(x)

  B   <- data.frame(matrix(par[(aik + 1):dTheta], p, k + 1))
  rownames(B)  <- colnames(x)
  colnames(B)  <- paste0("b", seq(0, k, 1))

  Ai     <- A |> dplyr::mutate(rowId = dplyr::row_number())
  Bi     <- as.matrix(B[, 1:s])
  listA  <- list()
  listB  <- list()
  indica <- seq_len(nrow(x))
  xtemp  <- x |> `rownames<-`(seq_len(nrow(x)))

  # --- Bootstrap loop ---
  for (i in seq_len(resamples)) {
    sample_ind  <- as.matrix(sample(indica, replace = TRUE))
    xb  <- xtemp[sample_ind, ]
    nb  <- nrow(xb); pb <- ncol(xb)
    aikb <- nb * k; bjkb <- pb * (k + 1)

    ran1        <- as.matrix(runif(aikb))
    sample_ind2 <- rbind(sample_ind, sample_ind + nb)
    par1        <- ran1[sample_ind2, ]
    par2        <- runif(bjkb)
    param       <- c(par1, par2)

    res.b <- run_optimr(param, xb, use_type)

    Bb   <- matrix(res.b$par[(aikb + 1):(aikb + bjkb)], pb, k + 1)
    Ab   <- matrix(res.b$par[1:aikb], nb, k)
    Ab_j <- data.frame(Ab, rowId = sample_ind)
    Ab_j <- dplyr::distinct(Ab_j, rowId, .keep_all = TRUE)

    if (sup) {
      Asup  <- dplyr::anti_join(Ai, Ab_j, by = "rowId")
      n_sup <- nrow(Asup)
      pars  <- runif(n_sup * k)
      xsup  <- xtemp[Asup$rowId, ]

      res.sup <- suppressWarnings(
        optimr(par = pars, fn = Indsup.BIN, gr = Indsup.GradBIN,
               xs = xsup, B = Bb, k = k, lambda = L, method = method)
      )

      Asup <- data.frame(matrix(res.sup$par, n_sup, k),
                         rowId = Asup[, ncol(Asup)])
      Abr  <- rbind(Ab_j, Asup) |> dplyr::arrange(rowId)
      Ab_m <- as.matrix(Abr |> dplyr::select(-rowId))

      outA <- procOPA(as.matrix(A), Ab_m, scale = TRUE, reflect = TRUE)
      Arot <- outA$Bhat
      R    <- outA$R
    } else {
      Abr  <- Ab_j |> dplyr::arrange(rowId)
      Ai_b <- dplyr::inner_join(Ai, Abr, by = "rowId") |>
        dplyr::select(dplyr::all_of(paste0("Dim", seq_len(k))))
      Ab_m <- as.matrix(Abr |> dplyr::select(-rowId))

      outA <- procOPA(as.matrix(Ai_b), Ab_m, scale = TRUE, reflect = TRUE)
      Arot <- outA$Bhat
      R    <- outA$R
    }

    Brot <- Bb[, 2:s] %*% R

    Arot <- data.frame(Arot, ind = Abr$rowId)
    colnames(Arot) <- c(paste0("Dim", seq_len(k)), "ind")
    rownames(Arot) <- rownames(x)[Abr$rowId]

    Brot <- data.frame(Bb[, 1], Brot, param = seq_len(p))
    colnames(Brot) <- c(paste0("b", seq(0, k, 1)), "param")
    rownames(Brot) <- rownames(B)

    Arot$resample <- i; Brot$resample <- i
    listA[[i]] <- Arot
    listB[[i]] <- Brot
  }

  ResultA <- dplyr::bind_rows(listA)
  ResultB <- dplyr::bind_rows(listB)
  rows    <- rownames(x)
  cols    <- colnames(x)

  # aggregate then rename so downstream columns are Dim1..k and bb0..k
  EspA_agg <- aggregate(. ~ ind,   ResultA, mean)
  EspB_agg <- aggregate(. ~ param, ResultB, mean)

  # EspA_agg columns: ind, Dim1..k, resample  (aggregate keeps original names)
  # We only need the Dim columns for Ahat
  dim_cols <- paste0("Dim", seq_len(k))
  bb_cols  <- paste0("bb",  seq(0, k, 1))

  Ahat <- EspA_agg[, dim_cols, drop = FALSE]          # keeps Dim1..k
  rownames(Ahat) <- rownames(x)

  # Rename EspB columns for the ellipse join and for Bhat
  colnames(EspB_agg) <- c("param",
                           paste0("bb", seq(0, k, 1)),
                           "resampleb")
  Bhat <- EspB_agg[, bb_cols, drop = FALSE]
  rownames(Bhat) <- colnames(x)

  # --- Confidence ellipses ---
  if (ellipses) {
    # For ellipses we need a "Dimb" version of EspA for the join
    EspA_join <- EspA_agg
    colnames(EspA_join) <- c("ind",
                              paste0("Dimb", seq_len(k)),
                              "resampleb")

    CentBootA <- dplyr::left_join(ResultA, EspA_join, by = "ind") |>
      dplyr::mutate(Dim1c = .data[["Dim1"]] - .data[["Dimb1"]],
                    Dim2c = .data[["Dim2"]] - .data[["Dimb2"]]) |>
      dplyr::select("ind", "resample", "Dim1c", "Dim2c")

    n_ell    <- nrow(EspA_join)
    datalist <- list()
    q        <- conf

    for (a in seq_len(n_ell)) {
      temp  <- dplyr::filter(CentBootA, .data[["ind"]] == a) |>
        dplyr::select("Dim1c", "Dim2c")
      sol   <- svd(temp)
      kk    <- sqrt(rowSums(sol$u * sol$u))
      r     <- quantile(kk, q)
      th    <- seq(0, 2 * pi, length.out = 500)
      z     <- rbind(r * cos(th), r * sin(th))
      esp   <- dplyr::filter(EspA_join, .data[["ind"]] == a) |>
        dplyr::select(-"ind", -"resampleb")
      v     <- as.data.frame(
        as.matrix(rep(1, 500)) %*% as.matrix(esp) +
          r * t(z) %*% diag(sol$d) %*% t(sol$v))
      colnames(v) <- paste0("Dim", seq_len(k))
      v$ind <- a
      datalist[[a]] <- v
    }
    Ellip <- dplyr::bind_rows(datalist)
  }

  # --- Confusion matrix ---
  n_orig <- nrow(x)
  LogP   <- as.matrix(cbind(rep(1, n_orig), Ahat)) %*% t(Bhat)
  Prob   <- plogis(LogP)
  Pr     <- ifelse(Prob >= 0.5, 1, 0)

  PCC   <- ifelse((x == 1 & Pr == 1) | (x == 0 & Pr == 0), 1, 0)
  ones  <- apply(x, 2, sum)
  zeros <- n_orig - ones

  confusion <- data.frame(
    Sensitivity = round(100 * apply((Pr == 1) & (x == 1), 2, sum) / ones,  1),
    Specificity = round(100 * apply((Pr == 0) & (x == 0), 2, sum) / zeros, 1),
    Global      = round(100 * colSums(PCC) / nrow(PCC), 1)
  )

  if (method == "CG" & type == 1) method <- "CG: Fletcher--Reeves"
  if (method == "CG" & type == 2) method <- "CG: Polak--Ribiere"
  if (method == "CG" & type == 3) method <- "CG: Beale--Sorenson"
  if (method == "Rcgmin")          method <- "CG: Dai--Yuan"

  if (ellipses) {
    out <- list(Ahat = Ahat, Bhat = Bhat, BootA = ResultA, BootB = ResultB,
                pred = Pr, fit = confusion, rows = rows, cols = cols,
                method = method, Ellip = Ellip)
  } else {
    out <- list(Ahat = Ahat, Bhat = Bhat, BootA = ResultA, BootB = ResultB,
                pred = Pr, fit = confusion, rows = rows, cols = cols,
                method = method)
  }

  if (plot & ncol(Ahat) > 1) {
    print(plotBLB(x = out, ellipses = ellipses, endsegm = 0.95,
                  col.ind = col.ind))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

# Internal: redirect ALL console output (stdout + messages) to a temp file.
.silently <- function(expr) {
  tmp <- tempfile()
  con <- file(tmp, open = "wt")
  sink(con, type = "output",  append = FALSE)
  sink(con, type = "message", append = TRUE)
  tryCatch(
    { result <- force(expr); result },
    error = function(e) { stop(e) },
    finally = {
      try(sink(type = "output"),  silent = TRUE)
      try(sink(type = "message"), silent = TRUE)
      try(close(con),             silent = TRUE)
      try(unlink(tmp),            silent = TRUE)
    }
  )
}

# Internal: run ONE bootstrap resample. Identical logic to the sequential loop.
# Returns list(Arot, Brot) or NULL on failure.
# Self-contained so it can run on a parallel worker.
.boot_one <- function(seed_i, xtemp, indica, n_orig, p_orig, k, s, L,
                      method, type, use_type, A_mat, Ai_df, sup,
                      rn_x, rn_B) {
  set.seed(seed_i)
  tryCatch({
    n <- n_orig; p <- p_orig
    aik <- n * k; bjk <- p * (k + 1)

    sample_ind <- as.matrix(sample(indica, replace = TRUE))
    xb  <- xtemp[sample_ind, ]
    n   <- nrow(xb); p <- ncol(xb)
    aik <- n * k; bjk <- p * (k + 1)

    ran1        <- as.matrix(runif(aik))
    sample_ind2 <- rbind(sample_ind, sample_ind + n)
    par1        <- ran1[sample_ind2, ]
    par2        <- runif(bjk)
    param       <- c(par1, par2)

    if (use_type) {
      res.b <- .silently(suppressWarnings(
        optimr(par = param, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
               xt = xb, k = k, lambda = L, method = method,
               control = list(type = type))
      ))
    } else {
      res.b <- .silently(suppressWarnings(
        optimr(par = param, fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
               xt = xb, k = k, lambda = L, method = method)
      ))
    }
    for (f in c(file.path(tempdir(), "badhess.txt"),
                file.path(getwd(),   "badhess.txt")))
      if (file.exists(f)) unlink(f, force = TRUE)

    Bb   <- matrix(res.b$par[(aik + 1):(aik + bjk)], p, k + 1)
    Ab   <- matrix(res.b$par[1:aik], n, k)
    Ab_j <- data.frame(Ab, rowId = sample_ind)
    Ab_j <- Ab_j[!duplicated(Ab_j$rowId), ]

    if (sup) {
      Asup  <- Ai_df[!Ai_df$rowId %in% Ab_j$rowId, ]
      n_sup <- nrow(Asup)
      pars  <- runif(dim(Asup)[1] * (dim(Asup)[2] - 1))
      xsup  <- xtemp[Asup$rowId, ]

      res.sup <- .silently(suppressWarnings(
        optimr(par = pars, fn = Indsup.BIN, gr = Indsup.GradBIN,
               xs = xsup, B = Bb, k = k, lambda = L, method = method)
      ))
      for (f in c(file.path(tempdir(), "badhess.txt"),
                  file.path(getwd(),   "badhess.txt")))
        if (file.exists(f)) unlink(f, force = TRUE)

      Asup <- data.frame(matrix(res.sup$par, n_sup, k),
                         rowId = Asup[, ncol(Asup)])
      Abr  <- rbind(Ab_j, Asup) |> dplyr::arrange(rowId)
      Ab   <- Abr |> dplyr::select(-rowId)
      outA <- procOPA(A_mat, as.matrix(Ab), scale = TRUE, reflect = TRUE)

    } else {
      Abr  <- Ab_j |> dplyr::arrange(rowId)
      Ai_b <- dplyr::inner_join(
                data.frame(A_mat, rowId = seq_len(nrow(A_mat))),
                Abr, by = "rowId") |>
        dplyr::select(dplyr::all_of(paste0("Dim", seq_len(k))))
      Ab   <- Abr |> dplyr::select(-rowId)
      outA <- procOPA(as.matrix(Ai_b), as.matrix(Ab), scale = TRUE, reflect = TRUE)
    }

    Arot <- data.frame(outA$Bhat, ind = Abr$rowId)
    colnames(Arot) <- c(paste0("Dim", seq_len(k)), "ind")
    rownames(Arot) <- rn_x[Abr$rowId]

    Brot <- data.frame(Bb[, 1], Bb[, 2:s] %*% outA$R,
                       param = seq_len(p_orig))
    colnames(Brot) <- c(paste0("b", seq(0, k, 1)), "param")
    rownames(Brot) <- rn_B

    list(Arot = Arot, Brot = Brot)
  }, error = function(e) NULL)
}



#' @importFrom optimx optimr
#' @importFrom shapes procOPA
#' @importFrom stats runif
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom parallel clusterEvalQ parLapplyLB
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
#' Set \code{ncores > 1} to distribute the bootstrap resamples across multiple
#' CPU cores. On Windows a PSOCK cluster is used; on Unix/Mac a FORK cluster
#' is used instead (lower overhead, faster). Use \code{ncores = NULL} to
#' automatically use all available cores minus one.
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
#' @param ncores Number of CPU cores for parallel execution. Default \code{1}
#'   (sequential with progress bar). Set to \code{NULL} to use all available
#'   cores minus one, or specify an integer \eqn{\geq 2}.
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
#' \dontrun{
#' data("Methylation")
#' set.seed(02052020)
#' out_sup  <- bootBLB(x = Methylation, ellipses = FALSE)
#' out_boot <- bootBLB(x = Methylation, sup = FALSE, ellipses = TRUE)
#'
#' # Parallel execution: distribute resamples across multiple CPU cores.
#' # Use ncores = NULL to automatically detect available cores minus one,
#' # or specify a number explicitly. Parallel mode does not show a progress
#' # bar. Results are numerically identical to the sequential mode.
#'
#' # Check available cores on your machine
#' parallel::detectCores()
#'
#' # Use 4 cores (approx. 4x faster than sequential)
#' data("Methylation")
#' set.seed(02052020)
#' out_par <- bootBLB(x = Methylation, method = "CG", type = 1,
#'                    ellipses = TRUE, ncores = 4)
#'
#' # Use all available cores minus one
#' out_par <- bootBLB(x = Methylation, method = "CG", type = 1,
#'                    ellipses = TRUE, ncores = NULL)
#' }

bootBLB <- function(x, k = 2, L = 0, method = "CG", type = 1, plot = TRUE,
                    sup = TRUE, ellipses = FALSE, maxit = NULL,
                    resamples = 100, conf = 0.9, col.ind = NULL,
                    ncores = 1) {

  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required. Please install it.", call. = FALSE)

  n <- nrow(x); p <- ncol(x)
  aik <- n * k; bjk <- p * (k + 1)
  dTheta <- aik + bjk; s <- k + 1
  use_type <- (method == "CG")

  # --- Initial fit on full data ---
  message("Fitting initial logistic biplot...")

  if (use_type) {
    res <- .silently(suppressWarnings(
      optimr(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
             xt = x, k = k, lambda = L, method = method,
             control = list(type = type))
    ))
  } else {
    res <- .silently(suppressWarnings(
      optimr(par = runif(dTheta), fn = J.BipLog.BIN, gr = Grad.BipLog.BIN,
             xt = x, k = k, lambda = L, method = method)
    ))
  }
  for (f in c(file.path(tempdir(), "badhess.txt"),
              file.path(getwd(),   "badhess.txt")))
    if (file.exists(f)) unlink(f, force = TRUE)

  par <- res$par
  A   <- data.frame(matrix(res$par[1:aik], n, k))
  colnames(A)  <- paste0("Dim", seq_len(k))
  rownames(A)  <- rownames(x)

  B   <- data.frame(matrix(par[(aik + 1):dTheta], p, k + 1))
  rownames(B)  <- colnames(x)
  colnames(B)  <- paste0("b", seq(0, k, 1))

  Ai     <- A |> dplyr::mutate(rowId = dplyr::row_number())
  indica <- seq_len(n)
  xtemp  <- x |> `rownames<-`(seq_len(n))
  A_mat  <- as.matrix(A)

  # Pre-generate seeds for reproducibility across sequential and parallel modes
  seeds <- sample.int(.Machine$integer.max, resamples)

  # Resolve number of cores
  n_cores <- if (is.null(ncores)) max(1L, detectCores() - 1L) else
             max(1L, as.integer(ncores))

  # -----------------------------------------------------------------------
  # Bootstrap: parallel (n_cores > 1) or sequential with progress bar
  # -----------------------------------------------------------------------
  if (n_cores > 1) {

    # PSOCK cluster works on all platforms including Windows
    cl <- makeCluster(n_cores, type = "PSOCK")
    on.exit(stopCluster(cl), add = TRUE)

    # Export package functions needed on workers
    clusterExport(cl, varlist = c(
      ".silently", ".boot_one",
      "J.BipLog.BIN", "Grad.BipLog.BIN",
      "Indsup.BIN",   "Indsup.GradBIN"
    ), envir = environment())

    # Load required packages on each worker
    clusterEvalQ(cl, {
      library(optimx)
      library(shapes)
      library(dplyr)
    })

    # Export all data/parameters needed by .boot_one
    clusterExport(cl, varlist = c(
      "xtemp", "indica", "n", "p", "k", "s", "L",
      "method", "type", "use_type",
      "A_mat", "Ai", "sup"
    ), envir = environment())

    message("Running ", resamples, " bootstrap resamples on ",
            n_cores, " cores...")
    results <- parLapplyLB(cl, seeds, function(seed_i) {
      .boot_one(seed_i  = seed_i,
                xtemp   = xtemp,  indica  = indica,
                n_orig  = n,      p_orig  = p,
                k       = k,      s       = s,
                L       = L,      method  = method,
                type    = type,   use_type = use_type,
                A_mat   = A_mat,  Ai_df   = Ai,
                sup     = sup,
                rn_x    = rownames(x),
                rn_B    = rownames(B))
    })

  } else {

    message("Running bootstrap resamples...")
    pb      <- txtProgressBar(min = 0, max = resamples, style = 3, width = 50)
    results <- vector("list", resamples)

    for (i in seq_len(resamples)) {
      results[[i]] <- .boot_one(
        seed_i   = seeds[i],
        xtemp    = xtemp,  indica   = indica,
        n_orig   = n,      p_orig   = p,
        k        = k,      s        = s,
        L        = L,      method   = method,
        type     = type,   use_type = use_type,
        A_mat    = A_mat,  Ai_df    = Ai,
        sup      = sup,
        rn_x     = rownames(x),
        rn_B     = rownames(B))
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  # Attach resample index; drop any NULLs (failed resamples)
  ok <- !vapply(results, is.null, logical(1))
  results <- results[ok]
  if (length(results) == 0)
    stop("All bootstrap resamples failed. Try a different method or seed.")

  listA <- lapply(seq_along(results), function(i) {
    r <- results[[i]]$Arot; r$resample <- i; r })
  listB <- lapply(seq_along(results), function(i) {
    r <- results[[i]]$Brot; r$resample <- i; r })

  ResultA <- dplyr::bind_rows(listA)
  ResultB <- dplyr::bind_rows(listB)
  rows    <- rownames(x)
  cols    <- colnames(x)

  EspA_agg <- aggregate(. ~ ind,   ResultA, mean)
  EspB_agg <- aggregate(. ~ param, ResultB, mean)

  dim_cols <- paste0("Dim", seq_len(k))
  bb_cols  <- paste0("bb",  seq(0, k, 1))

  Ahat <- EspA_agg[, dim_cols, drop = FALSE]
  rownames(Ahat) <- rownames(x)

  colnames(EspB_agg) <- c("param", paste0("bb", seq(0, k, 1)), "resampleb")
  Bhat <- EspB_agg[, bb_cols, drop = FALSE]
  rownames(Bhat) <- colnames(x)

  # --- Confidence ellipses ---
  if (ellipses) {
    EspA_join <- EspA_agg
    colnames(EspA_join) <- c("ind", paste0("Dimb", seq_len(k)), "resampleb")

    CentBootA <- dplyr::left_join(ResultA, EspA_join, by = "ind") |>
      dplyr::mutate(Dim1c = .data[["Dim1"]] - .data[["Dimb1"]],
                    Dim2c = .data[["Dim2"]] - .data[["Dimb2"]]) |>
      dplyr::select("ind", "resample", "Dim1c", "Dim2c")

    n_ell    <- nrow(EspA_join)
    datalist <- list()
    q        <- conf

    for (a in seq_len(n_ell)) {
      temp <- dplyr::filter(CentBootA, .data[["ind"]] == a) |>
        dplyr::select("Dim1c", "Dim2c")

      if (nrow(temp) < 2) next

      sol   <- svd(temp)
      d_pad <- c(sol$d, rep(0, k - length(sol$d)))
      kk    <- sqrt(rowSums(sol$u * sol$u))
      r     <- quantile(kk, q)
      th    <- seq(0, 2 * pi, length.out = 500)
      z     <- rbind(r * cos(th), r * sin(th))
      esp   <- dplyr::filter(EspA_join, .data[["ind"]] == a) |>
        dplyr::select(-"ind", -"resampleb")
      v     <- as.data.frame(
        as.matrix(rep(1, 500)) %*% as.matrix(esp) +
          r * t(z) %*% diag(d_pad) %*% t(sol$v))
      colnames(v) <- paste0("Dim", seq_len(k))
      v$ind <- a
      datalist[[a]] <- v
    }
    Ellip <- dplyr::bind_rows(datalist[!vapply(datalist, is.null, logical(1))])
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

  method_label <- method
  if (method == "CG" & type == 1) method_label <- "CG: Fletcher--Reeves"
  if (method == "CG" & type == 2) method_label <- "CG: Polak--Ribiere"
  if (method == "CG" & type == 3) method_label <- "CG: Beale--Sorenson"
  if (method == "Rcgmin")          method_label <- "CG: Dai--Yuan"
  method <- method_label

  if (ellipses && exists("Ellip") && nrow(Ellip) > 0) {
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
                  col.ind  = col.ind,
                  titles   = "Bootstrap Binary Logistic Biplot",
                  subtitle = paste0("Estimation with ", method)))
  }

  class(out) <- c("BiplotML", "list")
  return(out)
}

#' @title
#' Plot a Binary Logistic Biplot using a BiplotML object
#' @description
#' Plot the bootstrap binary logistic biplot and draw confidence ellipses on the individuals of an object BiplotML.
#' @return
#' Returns the Biplot of the individuals and variables.
#' @details
#' If draw = "ind", then the biplot is plotted only for individuals and if draw = "var" then is plotted only for the variables.
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x Object class BiplotML.
#' @param dim Dimensions plot. By default \code{Dim1} and \code{Dim2}.
#' @param col.ind Color for the individuals.
#' @param col.var Color for the variables.
#' @param label.ind By default the row points are not labelled.
#' @param draw The graph to draw ("ind" for the individuals, "var" for the variables and "biplot" for the row and columns coordinates in the same graph)
#' @param titles Title for the Biplot
#' @param ellipses If \code{ellipses=TRUE}, draw confidence ellipses around the rows.
#' @param endsegm Represents where the segment of a variable ends on the logit probability scale. By default \code{endsegm=0.75}
#' @param repel  Repel overlapping text labels.
#' @param xylim vector specifying the minimum and maximum of the x-axis and y-axis. For example, you can use xylim=c(-10, 10).
#' @references
#' Meulman, J. J., & Heiser, W. J. (1983). The display of bootstrap solutions in multidimensional scaling. Murray Hill, NJ: Bell Laboratories. (Technical memorandum)
#'
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{bootBLB}}
#' @examples
#' \donttest{
#' data("Methylation")
#' set.seed(123456)
#' outBLB <- bootBLB(x = Methylation, sup = TRUE, plot=FALSE)
#' plotBLB(x = outBLB, titles = "Methylation Logistic Biplot", ellipses = FALSE)
#' plotBLB(x = outBLB, titles = "Methylation LogBiplot", endsegm = 0.95)
#' plotBLB(x = outBLB, label.ind = TRUE, titles = "Methylation LogBiplot")
#' }
#' @export

plotBLB <- function(x, dim=c(1, 2), col.ind = NULL, col.var="#0E185F",
                    label.ind = FALSE, draw = c("biplot","ind","var"), titles = NULL,
                    ellipses = FALSE, endsegm = 0.75, repel = FALSE,
                    xylim = NULL){

  EspA <- x$Ahat
  EspB <- x$Bhat

  k <- ncol(EspA)

  if(k < max(dim)){
    stop(paste("The maximum number of dimensions in the logistic biplot fit is", k, "Review the dimensions requested for the biplot."))
  }

  grap <- match.arg(draw[1], c("ind","var","biplot"))

  #######.... Markers row
  dd1 <- paste0("bb", dim[1])
  dd2 <- paste0("bb", dim[2])
  EspB$x.50 = (-EspB[["bb0"]]*EspB[[dd1]]) / rowSums(EspB[,c(dd1, dd2)]*EspB[,c(dd1, dd2)])
  EspB$y.50 = (-EspB[["bb0"]]*EspB[[dd2]]) / rowSums(EspB[,c(dd1, dd2)]*EspB[,c(dd1, dd2)])

  EspB$x.75 = (log(endsegm/(1-endsegm))-EspB[["bb0"]])*EspB[[dd1]] / rowSums(EspB[,c(dd1, dd2)]*EspB[,c(dd1, dd2)])
  EspB$y.75 = (log(endsegm/(1-endsegm))-EspB[["bb0"]])*EspB[[dd2]] / rowSums(EspB[,c(dd1, dd2)]*EspB[,c(dd1, dd2)])

  colnames(EspA) <- c(paste0("Dim", seq(1,k,1)))

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package \"tidyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package \"dplyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package \"ggrepel\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(ellipses){
    min.plot <- ifelse(is.null(xylim), floor(min(min(x$Ellip$Dimb1), min(x$Ellip$Dimb2))), xylim[1])
    max.plot <- ifelse(is.null(xylim), ceiling(max(max(x$Ellip$Dimb1), max(x$Ellip$Dimb2))), xylim[2])
  }else{
    min.plot <- ifelse(is.null(xylim), floor(min(min(EspA$Dim1), min(EspA$Dim2))), xylim[1])
    max.plot <- ifelse(is.null(xylim), ceiling(max(max(EspA$Dim1), max(EspA$Dim2))), xylim[2])
  }

  if(x$method == "CG") subt <- "Estimation with Conjugate Gradient"
  if(x$method != "CG") subt <- paste0("Estimation with ", x$method, " algorithm")



  if(grap=="ind"){

    if (is.null(titles)) titulo <- "Individuals plot"
    if (!is.null(titles)) titulo <- titles

    if(is.null(col.ind)){
      g <- ggplot2::ggplot() +
        ggplot2::geom_point(data=EspA, ggplot2::aes(EspA[,dim[1]], y=EspA[,dim[2]]), color='#E7B800', size = 2, shape=17)+
        ggplot2::geom_vline(xintercept = 0, color = "black", linetype = 2) +
        ggplot2::geom_hline(yintercept = 0, color = "black", linetype = 2)  +
        ggplot2::xlab(paste("Dimension", dim[1])) +
        ggplot2::ylab(paste("Dimension", dim[2])) +
        ggplot2::xlim(min.plot, max.plot) + ggplot2::ylim(min.plot, max.plot) +
        ggplot2::scale_color_brewer(palette="Set1") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(face = "bold"), legend.position = "none") +
        ggplot2::labs(title = titulo,
                      subtitle = subt,
                      caption = Sys.Date()) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle =  ggplot2::element_text(hjust = 0.5))
    }else{
      g <-  ggplot2::ggplot() +
        ggplot2::geom_point(data=EspA, ggplot2::aes(EspA[,dim[1]], y=EspA[,dim[2]], colour=col.ind), size = 2, shape=17)+
        ggplot2::geom_vline(xintercept = 0, color = "black", linetype = 2) +
        ggplot2::geom_hline(yintercept = 0, color = "black", linetype = 2)  +
        ggplot2::xlab(paste("Dimension", dim[1])) +
        ggplot2::ylab(paste("Dimension", dim[2])) +
        ggplot2::xlim(min.plot, max.plot) + ggplot2::ylim(min.plot, max.plot) +
        ggplot2::scale_color_brewer(palette="Set1") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(face = "bold")) +
        ggplot2::labs(title = titulo,
                      subtitle = subt,
                      caption = Sys.Date()) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle =  ggplot2::element_text(hjust = 0.5))
    }
    if(label.ind){
      if(repel){
        g <- g +
          ggrepel::geom_text_repel(data=EspA, ggplot2::aes(x=EspA[,dim[1]], y=EspA[,dim[2]], label = rownames(EspA)), size=3.5, segment.color = "grey50")
      }else{
        g <- g +
          ggplot2::geom_text(data=EspA, ggplot2::aes(x=EspA[,dim[1]], y=EspA[,dim[2]], label = rownames(EspA)), size=3.5)
      }

    }
  }

  min.plot2 <- ifelse(is.null(xylim), floor(min(min(EspB$x.50), min(EspB$x.75), min(EspB$y.50), min(EspB$y.75))), xylim[1])
  max.plot2 <- ifelse(is.null(xylim), ceiling(max(max(EspB$x.50), max(EspB$x.75), max(EspB$y.50), max(EspB$y.75))), xylim[2])
  limit2 <- max(abs(min.plot2), abs(max.plot2))

  if(grap == "var"){

    if (is.null(titles)) titulo <- "Variables plot"
    if (!is.null(titles)) titulo <- titles

    g <- ggplot2::ggplot() +
      ggplot2::geom_segment(data=EspB, ggplot2::aes(x=x.50, y=y.50, xend=x.75, yend=y.75),
                            arrow=grid::arrow(angle=25, type="closed",ends="last", length=grid::unit(0.3,"cm")),
                            colour = col.var, show.legend=FALSE) +
      ggplot2::xlab(paste("Dimension", dim[1])) +
      ggplot2::ylab(paste("Dimension", dim[2])) +
      ggplot2::xlim(min.plot2, max.plot2) + ggplot2::ylim(min.plot2, max.plot2) +
      ggplot2::geom_vline(xintercept = 0, color = "black", linetype = 2) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linetype = 2)  +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text = ggplot2::element_text(face = "bold"),
                     legend.position = "none") +
      ggplot2::labs(title = titulo,
                    subtitle = subt,
                    caption = Sys.Date()) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle =  ggplot2::element_text(hjust = 0.5))

    if(repel){
      g <- g +
        ggrepel::geom_text_repel(data=EspB, ggplot2::aes(x=x.75, y=y.75, label = rownames(EspB),
                                                         angle = (180/pi) * atan(y.75/x.75),
                                                         hjust = (1 - 2 * sign(x.75))/ 2),
                                 size=3, segment.color = "grey50", colour = "black", vjust=0)
    }else{
      g <- g +
        ggplot2::geom_text(data=EspB, ggplot2::aes(x=x.75, y=y.75, label = rownames(EspB),
                                                   angle = (180/pi) * atan(y.75/x.75),
                                                   hjust = (1 - 2 * sign(x.75))/ 2),
                           size=3, colour = "black", vjust=0)
    }
  }

  min.plot3 <-  min(min.plot, min.plot2)
  max.plot3 <-  max(max.plot, max.plot2)
  limit3 <- max(abs(min.plot3), abs(max.plot3))

  if(grap == "biplot"){
    if (is.null(titles)) titulo <- "Logistic Biplot"
    if (!is.null(titles)) titulo <- titles

    if(is.null(col.ind)){
      g <- ggplot2::ggplot() +
        ggplot2::geom_point(data=EspA, ggplot2::aes(x=EspA[,dim[1]], y=EspA[,dim[2]]), color='#E7B800', size = 2, shape=17)+
        ggplot2::geom_segment(data=EspB, ggplot2::aes(x=x.50, y=y.50, xend=x.75, yend=y.75), arrow=grid::arrow(angle=20,type="closed",ends="last", length=grid::unit(0.3,"cm")), colour = col.var, show.legend=FALSE) +
        ggplot2::geom_vline(xintercept = 0, color = "black", linetype = 2) +
        ggplot2::geom_hline(yintercept = 0, color = "black", linetype = 2) +
        ggplot2::xlab(paste("Dimension", dim[1])) +
        ggplot2::ylab(paste("Dimension", dim[2])) +
        ggplot2::xlim(min.plot3, max.plot3) + ggplot2::ylim(min.plot3, max.plot3) +
        ggplot2::scale_color_brewer(palette="Set1") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(face = "bold"),
                       legend.position = "none") +
        ggplot2::labs(title = titulo,
                      subtitle = subt,
                      caption = Sys.Date()) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle =  ggplot2::element_text(hjust = 0.5))
    }else{
      g <- ggplot2::ggplot() +
        ggplot2::geom_point(data=EspA, ggplot2::aes(x=EspA[,dim[1]], y=EspA[,dim[2]], colour=col.ind), size = 2, shape=17)+
        ggplot2::geom_segment(data=EspB, ggplot2::aes(x=x.50, y=y.50, xend=x.75, yend=y.75), arrow=grid::arrow(angle=20,type="closed",ends="last", length=grid::unit(0.3,"cm")), colour = col.var, show.legend=FALSE) +
        ggplot2::geom_vline(xintercept = 0, color = "black", linetype = 2) +
        ggplot2::geom_hline(yintercept = 0, color = "black", linetype = 2) +
        ggplot2::xlab(paste("Dimension", dim[1])) +
        ggplot2::ylab(paste("Dimension", dim[2])) +
        ggplot2::xlim(min.plot3, max.plot3) + ggplot2::ylim(min.plot3, max.plot3) +
        ggplot2::scale_color_brewer(palette="Set1") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(face = "bold")) +
        ggplot2::labs(title = titulo,
                      subtitle = subt,
                      caption = Sys.Date()) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle =  ggplot2::element_text(hjust = 0.5))
    }

    if(repel){
      g <- g +
        ggrepel::geom_text_repel(data=EspB, ggplot2::aes(x=x.75, y=y.75, label = rownames(EspB),
                                                         angle = (180/pi) * atan(y.75/x.75), hjust = (1 - 2 * sign(x.75))/ 2),
                                 size=3, segment.color = "grey50", vjust=0)
    }else{
      g <- g +
        ggplot2::geom_text(data=EspB, ggplot2::aes(x=x.75, y=y.75, label = rownames(EspB),
                                                   angle = (180/pi) * atan(y.75/x.75), hjust = (1 - 2 * sign(x.75))/ 2),
                           size=3, vjust=0)
    }

    if(label.ind){
      g <- g +
        ggrepel::geom_text_repel(data=EspA, ggplot2::aes(x=Dim1, y=Dim2, label = rownames(EspA)), size=3.5, segment.color = "grey50")
    }
  }

  if(ellipses){
    g <- g +
      ggplot2::geom_point(data=x$Ellip, ggplot2::aes(x=Dimb1, y=Dimb2, group=ind), size=0.001, colour="lightgray", shape=20)
  }

  return(g + ggplot2::coord_fixed())
}

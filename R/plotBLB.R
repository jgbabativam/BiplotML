#' @import dplyr
#' @import ggplot2
#' @import ggrepel
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
#' @param label.ind By default the row points are labelled.
#' @param draw The graph to draw ("ind" for the individuals, "var" for the variables and "biplot" for the row and columns coordinates in the same graph)
#' @param titles Title for the Biplot
#' @param ellipses If \code{ellipses=TRUE}, draw confidence ellipses around the rows.
#' @param endsegm Represents where the segment of a variable ends on the logit probability scale. By default \code{endsegm=0.75}
#' @references
#' Meulman, J. J., & Heiser, W. J. (1983). The display of bootstrap solutions in multidimensional scaling. Murray Hill, NJ: Bell Laboratories. (Technical memorandum)
#' Vicente-Villardon, J.L. and Galindo, M. Purificacion (2006), \emph{Multiple Correspondence Analysis and related Methods. Chapter: Logistic Biplots}. Chapman-Hall
#' @seealso \code{\link{bootBLB}}
#' @examples
#' data("Methylation")
#' set.seed(123456)
#'
#' outBLB = bootBLB(x = Methylation, sup = TRUE, plot=FALSE)
#' plotBLB(x = outBLB, titles = "Methylation Logistic Biplot", ellipses = FALSE)
#' plotBLB(x = outBLB, titles = "Methylation LogBiplot", endsegm = 0.95)
#' plotBLB(x = outBLB, label.ind = TRUE, titles = "Methylation LogBiplot")
#' @export

plotBLB <- function(x, dim=c(1, 2), col.ind = "#E7B800", col.var="#00AFBB",
                     label.ind = FALSE, draw = c("biplot","ind","var"), titles = NULL,
                     ellipses = FALSE, endsegm = 0.75){

  EspA <- x$Ahat
  EspB <- x$Bhat

  k <- ncol(EspA)
  grap <- match.arg(draw[1], c("ind","var","biplot"))

  #######.... Markers row
  EspB$x.50 = (-EspB$bb0*EspB$bb1) / rowSums(EspB[,c("bb1", "bb2")]*EspB[,c("bb1", "bb2")])
  EspB$y.50 = (-EspB$bb0*EspB$bb2) / rowSums(EspB[,c("bb1", "bb2")]*EspB[,c("bb1", "bb2")])

  EspB$x.75 = (log(endsegm/(1-endsegm))-EspB$bb0)*EspB$bb1 / rowSums(EspB[,c("bb1", "bb2")]*EspB[,c("bb1", "bb2")])
  EspB$y.75 = (log(endsegm/(1-endsegm))-EspB$bb0)*EspB$bb2 / rowSums(EspB[,c("bb1", "bb2")]*EspB[,c("bb1", "bb2")])

  colnames(EspA) <- c(paste0("Dim", seq(1,k,1)))

  if(ellipses){
    min.plot <- floor(min(min(x$Ellip$Dimb1), min(x$Ellip$Dimb2)))
    max.plot <- ceiling(max(max(x$Ellip$Dimb1), max(x$Ellip$Dimb2)))
  }else{
    min.plot <- floor(min(min(EspA$Dim1), min(EspA$Dim2)))
    max.plot <- ceiling(max(max(EspA$Dim1), max(EspA$Dim2)))
  }

  if(x$method == "CG") subt <- "Estimation with Conjugate Gradient"
  if(x$method != "CG") subt <- paste0("Estimation with ", x$method, " algorithm")

  if(grap=="ind"){

    if (is.null(titles)) titulo <- "Individuals plot"
    if (!is.null(titles)) titulo <- titles

    g <-  ggplot() +
      geom_point(data=EspA, aes(x=Dim1, y=Dim2), colour=col.ind, size = 2, shape=17)+
      geom_vline(xintercept = 0, color = "black", linetype = 2) +
      geom_hline(yintercept = 0, color = "black", linetype = 2)  +
      xlab("Dimension 1") + ylab("Dimension 2") +
      geom_text_repel(data=EspA, aes(x=Dim1, y=Dim2, label = rownames(EspA)), size=3.5, segment.color = "grey50") +
      xlim(min.plot, max.plot) + ylim(-min.plot, max.plot) +
      theme_bw() +
      theme(axis.text = element_text(face = "bold"), legend.position = "none") +
      labs(title = titulo,
           subtitle = subt,
           caption = Sys.Date()) +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle =  element_text(hjust = 0.5))
  }

  min.plot2 <- floor(min(min(EspB$x.50), min(EspB$x.75), min(EspB$y.50), min(EspB$y.75)))
  max.plot2 <- ceiling(max(max(EspB$x.50), max(EspB$x.75), max(EspB$y.50), max(EspB$y.75)))
  limit2 <- max(abs(min.plot2), abs(max.plot2))

  if(grap == "var"){

    if (is.null(titles)) titulo <- "Variables plot"
    if (!is.null(titles)) titulo <- titles

    g <- ggplot() +
      geom_segment(data=EspB, aes(x=x.50, y=y.50, xend=x.75, yend=y.75), arrow=arrow(angle=20,type="closed",ends="last", length=unit(0.3,"cm")), colour = col.var, show.legend=FALSE) +
      geom_text_repel(data=EspB, aes(x=x.75, y=y.75, label = rownames(EspB),
                                     angle = (180/pi) * atan(y.75/x.75), hjust = (1 - 2 * sign(x.75))/ 2),
                      size=3, segment.color = "grey50", vjust=0) +
      xlab("Dimension 1") + ylab("Dimension 2") +
      xlim(min.plot2, max.plot2) + ylim(min.plot2, max.plot2) +
      geom_vline(xintercept = 0, color = "black", linetype = 2) +
      geom_hline(yintercept = 0, color = "black", linetype = 2)  +
      theme_bw() +
      theme(axis.text = element_text(face = "bold"), legend.position = "none") +
      labs(title = titulo,
           subtitle = subt,
           caption = Sys.Date()) +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle =  element_text(hjust = 0.5))

  }

  min.plot3 <-  min(min.plot, min.plot2)
  max.plot3 <-  max(max.plot, max.plot2)
  limit3 <- max(abs(min.plot3), abs(max.plot3))

  if(grap == "biplot"){
    if (is.null(titles)) titulo <- "Bootstrap Binary Logistic Biplot"
    if (!is.null(titles)) titulo <- titles

    g <- ggplot() +
      geom_point(data=EspA, aes(x=Dim1, y=Dim2), colour=col.ind, size = 2, shape=17)+
      geom_segment(data=EspB, aes(x=x.50, y=y.50, xend=x.75, yend=y.75), arrow=arrow(angle=20,type="closed",ends="last", length=unit(0.3,"cm")), colour = col.var, show.legend=FALSE) +
      geom_text_repel(data=EspB, aes(x=x.75, y=y.75, label = rownames(EspB),
                                     angle = (180/pi) * atan(y.75/x.75), hjust = (1 - 2 * sign(x.75))/ 2),
                      size=3, segment.color = "grey50", vjust=0) +
      geom_vline(xintercept = 0, color = "black", linetype = 2) +
      geom_hline(yintercept = 0, color = "black", linetype = 2) +
      xlab("Dimension 1") + ylab("Dimension 2") +
      xlim(min.plot3, max.plot3) + ylim(min.plot3, max.plot3) +
      theme_bw() +
      theme(axis.text = element_text(face = "bold"), legend.position = "none") +
      labs(title = titulo,
           subtitle = subt,
           caption = Sys.Date()) +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle =  element_text(hjust = 0.5))

    if(label.ind){
      g <- g +
        geom_text_repel(data=EspA, aes(x=Dim1, y=Dim2, label = rownames(EspA)), size=3.5, segment.color = "grey50")
    }
  }

  if(ellipses){
    g <- g +
      geom_point(data=x$Ellip, aes(x=Dimb1, y=Dimb2, group=ind), size=0.001, colour="lightgray", shape=20)
  }
  return(g)
}

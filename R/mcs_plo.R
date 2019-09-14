#' Create MCS summary plot
#'
#' @param mcsres output from \code{\link{mcs_fun}}
#' @param alpha numeric value for transparency of background fill
#' @param xmax numeric value for upper limit of x-axis, defaults to 95th percentile of all linkage values
#'
#' @export
#'
#' @import ggplot2
mcs_plo <- function(mcsres, alpha = 0.4, xmax = NULL){

  # input
  toplo <- mcsres 
  if(is.null(xmax))
    xmax <- quantile(toplo$sitsedlnk, 0.95)[[1]]

  # plot
  p <- toplo %>%
    ggplot() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) +
    annotate("rect", xmin = 0, xmax = xmax, ymin = 0,    ymax = 0.25, alpha = alpha, fill = "red") +
    annotate("rect", xmin = 0, xmax = xmax, ymin = 0.25, ymax = 0.5,  alpha = alpha, fill = "orange") +
    annotate("rect", xmin = 0, xmax = xmax, ymin = 0.5,  ymax = 0.75, alpha = alpha, fill = "blue") +
    annotate("rect", xmin = 0, xmax = xmax, ymin = 0.75, ymax = 1,    alpha = alpha, fill = "darkgreen") +
    stat_ecdf(aes(sitsedlnk, color = contam),
              geom = "line", size = 1.5) +
    geom_vline(xintercept = 0.5, lty = "dashed", color = "blue") +
    annotate("text", x = xmax - (0.1 * xmax), y = .875, label = "Very Low", size = 6, colour = 'darkgreen') +
    annotate("text", x = xmax - (0.1 * xmax), y = .615, label = "Low", size = 6, colour = 'blue') +
    annotate("text", x = xmax - (0.1 * xmax), y = .375, label = "Moderate", size = 6, colour = 'orange') +
    annotate("text", x = xmax - (0.1 * xmax), y = .125, label = "High", size = 6, colour = 'red') +
    labs(x = 'Site Linkage Factor',
         y = 'Cumulative Proportion',
         title = 'Site Linkage',
         fill = '') +
    scale_colour_manual(values = c('darkblue', 'darkred', 'darkgreen', 'purple')) +
    theme_bw(base_family = 'serif', base_size = 16) +
    theme(legend.position="top", legend.title = element_blank(), legend.text = element_text(size = 14)) + 
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1))

  return(p)
  
}
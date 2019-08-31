#' Create MCS summary plot
#'
#' @param mcsres output from \code{\link{mcs_fun}}
#'
#' @export
#'
#' @import ggplot2
mcs_plo <- function(mcsres){

  # input
  toplo <- mcsres 
  nFactor <- quantile(toplo$sitsedlnk, 0.95)[[1]]
  al <- 0.4
  
  # plot
  p <- toplo %>%
    ggplot() +
    scale_x_continuous(expand = c(0, 0), limits = c(0,nFactor)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    annotate("rect", xmin = 0, xmax = nFactor, ymin = 0,    ymax = 0.25, alpha = al, fill = "red") +
    annotate("rect", xmin = 0, xmax = nFactor, ymin = 0.25, ymax = 0.5,  alpha = al, fill = "orange") +
    annotate("rect", xmin = 0, xmax = nFactor, ymin = 0.5,  ymax = 0.75, alpha = al, fill = "blue") +
    annotate("rect", xmin = 0, xmax = nFactor, ymin = 0.75, ymax = 1,    alpha = al, fill = "darkgreen") +
    stat_ecdf(aes(sitsedlnk, color = contam),
              geom = "line", size = 1.5) +
    geom_vline(xintercept = 0.5, lty = "dashed", color = "blue") +
    annotate("text", x = nFactor-1, y = .875, label = "Very Low", size = 6, colour = 'darkgreen') +
    annotate("text", x = nFactor-1, y = .615, label = "Low", size = 6, colour = 'blue') +
    annotate("text", x = nFactor-1, y = .375, label = "Moderate", size = 6, colour = 'orange') +
    annotate("text", x = nFactor-1, y = .125, label = "High", size = 6, colour = 'red') +
    labs(x = 'Site Linkage Factor',
         y = 'Cumulative Proportion',
         title = 'Site Linkage',
         fill = '') +
    scale_colour_manual(values = c('darkblue', 'darkred', 'darkgreen', 'purple')) +
    theme_bw(base_family = 'serif', base_size = 16) +
    theme(legend.position="top", legend.title = element_blank(), legend.text = element_text(size = 14)) 

  return(p)
  
}
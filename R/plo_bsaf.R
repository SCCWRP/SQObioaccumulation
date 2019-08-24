#' plot bsaf and cbiota for contaminant
#'
#' @param bsaf bsaf output from \code{\link{bioaccum_batch}}
#' @param cbiota cbiota output from \code{\link{bioaccum_batch}} 
#' @param cntbsaf selected contaminant
#'
#' @export
#' 
#' @import ggplot2
#'
plo_bsaf <- function(bsaf, cbiota, cntbsaf){
  
  # combine for plotting
  toplo <- rescmb(bsaf, cbiota, cntbsaf)
  
  p <- ggplot(toplo, aes(x = species, y = val)) + 
    geom_bar(stat = 'identity') + 
    xlab(cntbsaf) + 
    facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y') + 
    theme_bw(base_size= 16, base_family = 'serif') + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_blank(),
      strip.background = element_blank(), 
      strip.placement = 'outside'
    )
  
  return(p)
  
}
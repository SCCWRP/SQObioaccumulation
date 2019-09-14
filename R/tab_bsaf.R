#' Table of bsaf and cbiota for contaminant
#'
#' @param bsaf bsaf output from \code{\link{bioaccum_batch}}
#' @param cbiota cbiota output from \code{\link{bioaccum_batch}} 
#' @param cntbsaf selected contaminant
#'
#' @export
#'
#' @import dplyr tidyr
#' @importFrom magrittr "%>%"
tab_bsaf <- function(bsaf, cbiota, cntbsaf = NULL){
  
  # combine for plotting
  totab <- rescmb(bsaf, cbiota, cntbsaf) %>% 
    mutate(val = round(val, 3)) %>% 
    spread(species, val) %>% 
    rename(Output = var) %>% 
    mutate(Output = factor(Output, levels = c('Tissue conc. (ng/g wet)', 'BSAF'))) %>% 
    arrange(Output)
  
  return(totab)
  
}
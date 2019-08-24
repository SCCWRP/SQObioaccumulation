#' Combine bsaf and cbiota for plot/table of contaminant
#'
#' @param bsaf bsaf output from \code{\link{bioaccum_batch}}
#' @param cbiota cbiota output from \code{\link{bioaccum_batch}} 
#' @param cntbsaf selected contaminant
#'
#' @export
#' 
#' @import dplyr tidyr
#' @importFrom magrittr "%>%"
rescmb <- function(bsaf, cbiota, cntbsaf){
  
  # check if selected contaminant in bsaf
  if(!cntbsaf %in% names(bsaf)[!names(bsaf) %in% 'species'])
    stop('contaminant not found, must be in bsaf input')  
  
  # bsaf
  bsaf <- bsaf %>% 
    select(species, !!cntbsaf) %>% 
    rename(
      BSAF = !!cntbsaf
    ) %>% 
    mutate(species = factor(species, levels = species))
  
  # cbiota
  cbiota <- cbiota %>% 
    select(species, !!cntbsaf) %>% 
    rename(
      `Tissue conc. (ng/g wet)` = !!cntbsaf
    ) %>% 
    mutate(species = factor(species, levels = species))
  
  out <- bsaf %>% 
    left_join(cbiota, by = 'species') %>% 
    gather('var', 'val', -species)
  
  return(out)
  
}
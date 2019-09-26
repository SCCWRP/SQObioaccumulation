#' Combine bsaf and cbiota for plot/table of contaminant
#'
#' @param bsaf bsaf output from \code{\link{bioaccum_batch}}
#' @param cbiota cbiota output from \code{\link{bioaccum_batch}} 
#' @param cntbsaf selected contaminant, optional
#'
#' @export
#' 
#' @import dplyr tidyr
#' @importFrom magrittr "%>%"
rescmb <- function(bsaf, cbiota, cntbsaf = NULL){
  
  # species and contaminant levels
  splev <- bsaf$species
  cnlev <- names(cbiota)[!names(cbiota) %in% 'species']

  # bsaf
  bsaf <- bsaf %>% 
    pivot_longer(-species, names_to = 'Contaminant', values_to = 'BSAF')
  
  # cbiota
  cbiota <- cbiota %>%
    pivot_longer(-species, names_to = 'Contaminant', values_to = 'Tissue conc. (ng/g wet)')
  
  # combine both
  out <- bsaf %>% 
    left_join(cbiota, by = c('species', 'Contaminant')) %>% 
    pivot_longer(c(BSAF, `Tissue conc. (ng/g wet)`), names_to = 'var', values_to = 'val') %>% 
    mutate(
      species = factor(species, levels = splev), 
      Contaminant = factor(Contaminant, levels = cnlev)
      ) %>% 
    arrange(Contaminant)
  
  # return results for one 
  if(!is.null(cntbsaf)){
     
    # check if selected contaminant in bsaf
    if(!cntbsaf %in% out$Contaminant)
      stop('contaminant not found, must be in bsaf input')
    
    out <- out %>% 
      filter(Contaminant %in% cntbsaf) %>% 
      select(-Contaminant)
    
  }
 
  return(out)
  
}
#' site use mcs function
#'
#' @param nsim number of simulations
#' @param constants constants inputs
#' @param mcsparms MCS parameter inputs
#'
#' @export
#'
suf_mcs_fun <- function(nsim, constants, mcsparms){
  
  # site area and length
  SA <- constants %>% 
    filter(Constant %in% 'SA') %>% 
    pull(Value)
  SL <- constants %>% 
    filter(Constant %in% 'SL') %>%
    pull(Value)
  
  # home range mean and sd for guild species
  hrvals <- mcsparms %>% 
    filter(grepl('^HR[0-9]', MCSvar)) %>% 
    rename(species = MCSvar) %>% 
    mutate(
      var = case_when(
        grepl('X$', species) ~ 'X', 
        grepl('SD$', species) ~ 'SD'
      ),
      species = gsub('^HR', 'indic', species),
      species = gsub('X$|SD$', '', species)
    ) %>% 
    spread(var, Value)
  
  # home range sims
  sufsims <- hrvals %>% 
    group_by(species) %>% 
    mutate(
      suf = purrr::map(list(species), function(...){
        
        # indic1, indic8, indic9
        if(grepl('1$|8$|9$', species))
          out <- genlognorm_fun(nsim, X, SD) %>% 
            mutate(
              sims = SL / sims,
              sims = ifelse(is.infinite(sims), 0, sims)
            )
        
        # indic2, indic3, indic4, indic5, indic7
        if(grepl('2$|3$|4$|5$|7$', species))
          out <- genlognorm_fun(nsim, X, SD) %>% 
            mutate(
              sims = SA / sims,
              sims = ifelse(is.infinite(sims), 0, sims)
            )
        
        # indic6
        if(grepl('6$', species)){
          out <- (SL * 1000) / pgamma(runif(nsim, 0, 1), shape = X, scale = SD) 
          simi <- seq(1:nsim)
          out <- tibble(i = simi, sims = out)
        }
        
        return(out)
        
      })
    ) %>% 
    dplyr::select(-SD, -X) %>% 
    unnest(suf) %>% 
    mutate(
      sims = pmin(1, sims)
    ) %>% 
    rename(suf = sims)
  
  return(sufsims)
  
}
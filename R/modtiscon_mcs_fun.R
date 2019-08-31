#' mcs function for modelled tissue concentration
#'
#' @param nsim number of simulations
#' @param meanse mean and se values from user input for each guild species and contaminant class
#' @param propseaf proportion of seafood diet
#'
#' @export
modtiscon_mcs_fun <- function(nsim, meanse, propseaf){
  
  # simulated tissue concentrations across guild species, all sims
  sims <- meanse %>%  
    group_by(MCSvar, contam) %>% 
    mutate(
      ests = purrr::pmap(list(nsim, X, SD), genlognorm_fun)
    ) %>% 
    dplyr::select(-SD, -X) %>% 
    unnest
  
  # weighted tissue concentrations across guilds for each contam, all sims
  out <- sims %>% 
    dplyr::group_by(i) %>% 
    nest() %>% 
    mutate(
      wgtave = purrr::map(data, function(x){
        
        sumprod <- x %>% 
          arrange(contam, MCSvar) %>% 
          mutate(
            sims = case_when(
              is.na(sims) ~ 0, 
              T ~ sims
            )
          ) %>% 
          group_by(contam) %>% 
          summarise(
            wgtave = sims %*% propseaf
          )
        
        return(sumprod)
        
      })
    ) %>% 
    dplyr::select(-data)
  
  return(out)
  
}
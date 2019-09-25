#' mcs function for sediment concentration
#'
#' @param nsim number of simulations
#' @param sedmeanse sediment mean and se values from user input for each contaminant class
#' @param propseaf proportion of seafood diet, output from formpropseaf
#' @param SUF site use factor
#' @param CVBAF bioaccumulation factor sd/mean from mcsparms
#' @param indic_sum indicator guild concentrations sums across contaminants
#'
#' @export
modsedcon_mcs_fun <- function(nsim, sedmeanse, propseaf, SUF, CVBAF, indic_sum){
  
  # simulated sediment concentrations, all contams
  sedsims <- sedmeanse %>%  
    group_by(contam) %>% 
    mutate(
      ests = purrr::pmap(list(nsim, X, SD), genlognorm_fun) 
    ) %>% 
    dplyr::select(-SD, -X) %>% 
    unnest(ests) %>% 
    dplyr::ungroup()
  
  # bioaccumulation sims
  biosims <- indic_sum %>% 
    tidyr::gather('contam', 'val', -species) %>% 
    filter(grepl('\\_calc$', contam)) %>% 
    mutate(
      contam = gsub('\\_calc$', '', contam)
    ) %>% 
    dplyr::group_by(species, contam) %>% 
    mutate(
      ests = purrr::pmap(list(nsim, val, CVBAF * val), genlognorm_fun)
    ) %>% 
    dplyr::select(-val) %>% 
    dplyr::ungroup()
  
  # combine sediment sims with biosims and SUF
  estcncsims <- biosims %>% 
    unnest(ests) %>% 
    full_join(sedsims, ., by = c('contam', 'i')) %>% 
    full_join(SUF, by = c('i', 'species')) %>% 
    mutate(
      estcnc = sims.x * suf * sims.y
    ) %>% 
    dplyr::select(-MCSvar, -sims.x, -sims.y, -suf)

  # weighted tissue concentrations across guilds for each contam, all sims
  out <- estcncsims %>% 
    dplyr::group_by(i, contam) %>% 
    arrange(i, contam, species) %>% 
    mutate(
      estcnc = case_when(
        sum(is.na(estcnc)) == 9 ~ NaN,
        is.na(estcnc) ~ 0, 
        T ~ estcnc
      )
    ) %>% 
    summarise(
      wgtave = estcnc %*% propseaf
    )
  
  return(out)
  
}

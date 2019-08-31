#' MCS function
#'
#' @param nsim number of MC sims
#' @param indic_sum output from indic_sum_fun
#' @param mcsparms MCS parameter inputs
#' @param constants constants inputs
#' @param propseaf numeric vector indicating proportion of human diet for each guild species, must sum to 1 and have length 9
#'
#' @export
mcs_fun <- function(nsim, indic_sum, mcsparms, constants, propseaf = c(0, 0.5, 0, 0, 0.5, 0, 0, 0, 0)){
  
  # sanity checks
  stopifnot(length(propseaf) == 9)
  stopifnot(sum(propseaf) == 1)
  
  ##
  # inputs 
  
  # CVBAF
  CVBAF <- mcsparms %>% 
    filter(MCSvar == 'CVBAF') %>% 
    pull
  
  # mean and se values from observed contaminants, from user inputs
  meanse <- mcsparms %>% 
    filter(grepl('^indic', MCSvar)) %>% 
    mutate(
      var = case_when(
        grepl('X$', MCSvar) ~ 'X', 
        grepl('SD$', MCSvar) ~ 'SD'
      ),
      contam = gsub('^indic[0-9]|X$|SD$', '', MCSvar), 
      MCSvar = gsub('(^indic[0-9]).*$', '\\1', MCSvar)
    ) %>% 
    spread(var, Value)
  
  # mean and se values for sediment contaminants, from user inputs
  sedmeanse <- mcsparms %>% 
    filter(grepl('^sed', MCSvar)) %>% 
    mutate(
      var = case_when(
        grepl('X$', MCSvar) ~ 'X', 
        grepl('SD$', MCSvar) ~ 'SD'
      ),
      contam = gsub('^sed|X$|SD$', '', MCSvar), 
      MCSvar = gsub('(^sed).*$', '\\1', MCSvar)
    ) %>% 
    spread(var, Value)
  
  ##
  # modeled tissue concentration for consumption risk, mcs
  # returns weighted concentrations across all sims
  modtiscon <- modtiscon_mcs_fun(nsim, meanse, propseaf)
  
  ##
  # site use function sims
  SUF <- suf_mcs_fun(nsim, constants, mcsparms)
  
  ##
  # rename columns in indic_sum
  indic_sum <- indic_sum %>% 
    rename(
      species = Guild,
      Chlordane_calc = `Chlordanes BSAF (calc)`,
      Dieldrin_calc = `Dieldrin BSAF (calc)`,
      DDT_calc = `DDTs BSAF (calc)`,
      PCB_calc = `PCBs BSAF (calc)`,
      Chlordane_conc = `Chlordanes Conc (ng/g)`,
      Dieldrin_conc = `Dieldrin Conc (ng/g)`,
      DDT_conc = `DDTs Conc (ng/g)`,
      PCB_conc = `PCBs Conc (ng/g)`
    )
  
  ##
  # modeled sediment contribution to tissue concentration, mcs
  # returns weighted concentrations across all sims
  modsedcon <- modsedcon_mcs_fun(nsim, sedmeanse, propseaf, SUF, CVBAF, indic_sum)
  
  ## 
  # combine modeled tissue and sediment concentrations to get site linkages
  out <- modtiscon %>% 
    full_join(modsedcon, by = 'i') %>% 
    tidyr::unnest() %>% 
    mutate(sitsedlnk = wgtave1 / wgtave) %>% 
    dplyr::select(-wgtave, -contam1, -wgtave1)
  
  return(out)
  
}
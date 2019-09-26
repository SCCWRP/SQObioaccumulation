#' MCS function
#'
#' @param nsim number of MC sims
#' @param indic_sum output from indic_sum_fun
#' @param mcsparms MCS parameter inputs
#' @param constants constants inputs
#'
#' @export
mcs_fun <- function(nsim, indic_sum, mcsparms, constants){

  ##
  # inputs 
  
  # propseaf
  propseaf <- mcsparms %>% 
    filter(grepl('^indic[0-9]propseaf', MCSvar)) %>% 
    mutate(
      Value = ifelse(is.na(Value), 0, Value)
    ) %>% 
    arrange(MCSvar) %>% 
    pull(Value)
  
  # sanity checks
  stopifnot(length(propseaf) == 9)
  stopifnot(sum(propseaf) == 1)
  
  # CVBAF
  CVBAF <- mcsparms %>% 
    filter(MCSvar == 'CVBAF') %>% 
    pull
  
  # mean and se values from observed contaminants, from user inputs
  meanse <- mcsparms %>% 
    filter(grepl('^indic', MCSvar) & !grepl('^indic[0-9]propseaf', MCSvar)) %>% 
    mutate(
      var = case_when(
        grepl('X$', MCSvar) ~ 'X', 
        grepl('SD$', MCSvar) ~ 'SD'
      ),
      contam = gsub('^indic[0-9]|X$|SD$', '', MCSvar), 
      MCSvar = gsub('(^indic[0-9]).*$', '\\1', MCSvar)
    ) %>% 
    pivot_wider(names_from = var, values_from = Value) %>% 
    mutate(
      SD = ifelse(!is.na(X) & is.na(SD), 0, SD)
    )
  
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
    pivot_wider(names_from = var, values_from = Value) %>% 
    mutate(
      SD = ifelse(!is.na(X) & is.na(SD), 0, SD)
    )
  
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
    full_join(modsedcon, by = c('i', 'contam')) %>% 
    mutate(sitsedlnk = wgtave.y / wgtave.x) %>% # sed / tis
    dplyr::select(-wgtave.x, -wgtave.y)
  
  return(out)
  
}
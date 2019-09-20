#' Format calculated contaminant concentrations
#'
#' @param contam input contaminants from formatted user inputs
#' @param constants input constants from formatted user inputs
#' 
#' @import dplyr tidyr
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
cntcalc <- function(contam, constants){
  
  # total organic carbon input (ocsed)
  ocsed <- constants %>% 
    filter(Constant %in% 'ocsed') %>% 
    pull(Value)
  
  # salinity input (Sal)
  Sal <- constants %>% 
    filter(Constant %in% 'salinity') %>% 
    pull(Value)
  
  # temperature
  Temp <- constants %>% 
    filter(Constant %in% 'T') %>% 
    pull(Value)
  
  # poc concentration in water
  xpoc <- constants %>% 
    filter(Constant %in% 'xpoc') %>% 
    pull(Value)
  
  # doc concentration in water
  xdoc <- constants %>% 
    filter(Constant %in% 'xdoc') %>% 
    pull(Value)
  
  # disequalibrium factor for poc partitioning
  dpoc <- constants %>% 
    filter(Constant %in% 'dpoc') %>% 
    pull(Value)
  
  # disequalibrium factor for doc partitioning
  ddoc <- constants %>% 
    filter(Constant %in% 'ddoc') %>% 
    pull(Value)
  
  # proportionality constant describing phase partitioning of poc
  alphapoc <- constants %>% 
    filter(Constant %in% 'alphapoc') %>% 
    pull(Value)
  
  # proportionality constant describing phase partitioning of doc
  alphadoc <- constants %>% 
    filter(Constant %in% 'alphadoc') %>% 
    pull(Value)
  
  # log kow temp corrected
  # log kow sal corrected (whic includes temp correction)
  # phi
  # calculated dissolved surface water concentration
  # free dissolved surface water concentration
  # calculated porewater concentration
  # log koc
  contam <- contam %>% 
    mutate(
      cs_ng.g = ifelse(is.na(cs_ng.g), 0, cs_ng.g),
      logkow_tempcor = round(log10(Kow), 2) - ((delt_uow / (log(10) * 0.0083145)) * ((1 / (273 + Temp)) - (1 / 298))),
      log_KowTS = log10(10 ^ logkow_tempcor* (10 ^ (0.0018 * LeBas_Molar_Volume * (0.5 * Sal / 35)))),
      phi = 1 / (1 + (xpoc * dpoc * alphapoc * 10 ^ log_KowTS) + (xdoc * ddoc * alphadoc * 10 ^ log_KowTS)),
      calc_cd_pg.l = ifelse(!is.na(cs_ng.g), 1e6 * (cs_ng.g / (ocsed * (0.35 * 10 ^ log_KowTS)) / 8), 0),
      free_cd_ng.ml = ifelse(is.na(cd_ng.g), calc_cd_pg.l, ifelse(cd_ng.g <= calc_cd_pg.l, cd_ng.g, calc_cd_pg.l)) / 1e6,
      calc_cp_ng.ml = ifelse(cp_ng.g > 0 & !is.na(cp_ng.g), cp_ng.g / 1e6, ((cs_ng.g / ocsed) / (0.35 * 10 ^ log_KowTS))), 
      log_koc = log10(0.35 * 10 ^ log_KowTS)
    )
  
  return(contam)
  
}
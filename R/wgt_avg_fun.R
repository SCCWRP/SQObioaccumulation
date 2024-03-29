#' Calculate weighted average observed tissue concentration (ng/g), from empirical data
#'
#' @param mcsparms input mcsparms data frame, observed average concentrations extracted
#'
#' @export
wgt_avg_fun <- function(mcsparms){
  
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

  # observed contaminants from user input, mean only
  contobs <- mcsparms %>% 
    filter(grepl('^indic.*X$', MCSvar)) %>% 
    mutate(
      contam = gsub('^indic[0-9](.*)X$', '\\1', MCSvar), 
      MCSvar = gsub('(^indic[0-9]).*$', '\\1', MCSvar), 
      Value = case_when(
        is.na(Value) ~ 0, 
        T ~ Value
      )
    ) %>% 
    arrange(contam, MCSvar)
  
  # weighted average observed tissue conc
  wgt_avg <- contobs %>% 
    group_by(contam) %>%
    summarise(
      wgt_obs = case_when(
        sum(Value) == 0 ~ NaN, 
        sum(Value) != 0 ~ Value %*% propseaf
      )
    )
  
  return(wgt_avg)
  
}
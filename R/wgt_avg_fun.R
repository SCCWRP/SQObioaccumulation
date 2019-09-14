#' Calculate weighted average observed tissue concentration (ng/g), from empirical data
#'
#' @param mcsparms input mcsparms data frame, observed average concentrations extracted
#' @param propseaf numeric vector indicating proportion of human diet for each guild species
#'
#' @export
wgt_avg_fun <- function(mcsparms, propseaf = c(0, 0.5, 0, 0, 0.5, 0, 0, 0, 0)){
  
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
      wgt_obs = Value %*% propseaf
    )
  
  return(wgt_avg)
  
}
#' Summarize MCS results, compare with observed
#'
#' @param mcsres MC results, output from \code{\link{mcs_fun}}
#'
#' @export
mcs_sum_fun <- function(mcsres){
  
  # get percentiles
  persitsed <- mcsres %>% 
    group_by(contam) %>% 
    nest %>% 
    mutate(
      percnt = purrr::map(data, function(x){
        
        prc <- quantile(x$sitsedlnk, c(0, .01, .05, .1, 0.25, .5, .75, 0.9, .95, .99, 1)) %>% 
          enframe 
        
        return(prc)
        
      })
    ) %>% 
    dplyr::select(-data) %>% 
    unnest %>% 
    mutate(name = factor(name, levels = c('0%', '1%', '5%', '10%', '25%', '50%', '75%', '90%', '95%', '99%', '100%'))) %>% 
    rename(
      percentile = name,
      Compound = contam
    ) %>% 
    spread(percentile, value)
  
  
  return(persitsed)
  
}

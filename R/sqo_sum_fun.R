#' SQO assessment summary
#'
#' @param wgtavg weighted average observed tissue concentrations from input, by contaminant category, output from \code{wgt_avg_fun}
#' @param mcsres output from \code{mcs_fun}
#' @param constants constants from user inputs and lookup table, only SCT is used (sediment linkage threshold)
#'
#' @details data used from lazy load include \code{\link{tischmthr}}, \code{\link{finalsiteassess}}
#' 
#' @export
sqo_sum_fun <- function(wgtavg, mcsres, constants){
  
  # category scores and labels, final labels
  levs <- c('1', '2', '3', '4', '5')
  labs <- c('Very Low', 'Low', 'Moderate', 'High', 'Very High')
  flabs <- c('Unimpacted', 'Likely Unimpacted', 'Possibly Impacted', 'Likely Impacted', 'Clearly Impacted')

  # sediment linkage threshold
  SCT <- constants %>% 
    filter(Constant %in% 'SCT') %>% 
    pull(Value)
  
  # thresholds in nested format for join
  tischmthr <- tischmthr %>% 
    group_by(contam) %>% 
    nest(.key = 'thr')
  
  # quartiles from MCSsum
  mcsres <- mcsres %>% 
    mcs_sum_fun %>% 
    gather('percentile', 'value', -Compound) %>% 
    rename(contam = Compound) %>% 
    filter(grepl('25|50|75', percentile))
  
  # combined data to get category outcomes
  cmb <- wgtavg %>% 
    full_join(mcsres, by = 'contam') %>% 
    spread(percentile, value) %>% 
    full_join(tischmthr, by = 'contam') 
  
  # get category outcomes
  # chmscr/chmlab - chemical exposure category score
  # lnkscr/lnklab - site linkage category score
  # sitscr/sitlab - final site assessment category score
  sums <- cmb %>% 
    mutate(
      wgt_est = wgt_obs * `50%`,
      chmscr = purrr::pmap(list(wgt_obs, thr), function(wgt_obs, thr){
        
        val <- thr %>% pull(val)
        scr <- 1 + findInterval(wgt_obs, val)
        
        return(scr)
        
      }),
      chmlab = factor(as.character(chmscr), levels = levs, labels = labs), 
      chmlab = as.character(chmlab)
    ) %>% 
    rowwise() %>% 
    mutate(
      lnkscr = ifelse(is.na(`25%`), NA, 4 - findInterval(SCT, c(`25%`, `50%`, `75%`))), 
      lnklab = factor(as.character(lnkscr), levels = levs, labels = labs), 
      lnklab = as.character(lnklab)
    ) %>% 
    unite('cmbscr', chmscr, lnkscr, sep = ', ', remove = F) %>% 
    mutate(
      sitscr = factor(cmbscr, levels = finalsiteassess[[1]], labels = finalsiteassess[[2]]), 
      sitscr = as.numeric(as.character(sitscr)), 
      sitlab = factor(sitscr, levels = levs, labels = flabs), 
      sitlab = as.character(sitlab)
    )
  
  # final formatting (no calcs)
  out <- sums %>% 
    select(-thr) %>% 
    unnest(chmscr) %>% 
    select(
      Compound = contam, 
      `Observed tissue conc. (ng/g)` = wgt_obs,
      `Chemical exposure category` = chmlab, 
      `Estimated tissue conc. (ng/g)` = wgt_est,
      `Site linkage 25%` = `25%`, 
      `Site linkage 50%` = `50%`, 
      `Site linkage 75%` = `75%`, 
      `Site linkage category` = lnklab, 
      `Site assessment category` = sitlab
    )
  
  return(out)
  
}
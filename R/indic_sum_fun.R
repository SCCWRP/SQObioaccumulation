#' Summary function for calculating bsaf and cbiota totals for each indicator species
#'
#' @param cbiota cbiota output from \code{\link{bioaccum_batch}} 
#' @param contamcalc contaminants fom inputs
#'
#' @import dplyr tidyr
#' 
#' @export
#' 
#' @importFrom magrittr "%>%"
indic_sum_fun <- function(cbiota, contamcalc){
  
  # cbiota long format
  cbiota_lng <- cbiota %>% 
    pivot_longer(-species, names_to = 'Chem', values_to = 'val')
  
  # contaminant inputs
  contams <- contamcalc %>% 
    dplyr::select(ChemGroup, Chem, cs_ng.g)
  
  # combine, add group, summarize by indic, est, chemgroup
  sumdat <- cbiota_lng %>% 
    filter(grepl('^indic', species)) %>% 
    left_join(contams, by = 'Chem') %>% 
    group_by(species, ChemGroup) %>% 
    summarise(
      calc = sum(val)/sum(cs_ng.g), 
      conc = sum(val)
    ) %>% 
    ungroup %>% 
    pivot_longer(c(calc, conc), names_to = 'var', values_to = 'val') %>% 
    unite('var', ChemGroup, var) %>% 
    pivot_wider(names_from = var, values_from = val) %>% 
    mutate_if(is.numeric, round, 2)
  
  # prettify
  out <- sumdat %>% 
    select(species, Chlordane_calc, Dieldrin_calc, DDT_calc, PCB_calc, Chlordane_conc, Dieldrin_conc, DDT_conc, PCB_conc) %>% 
    rename(
      Guild = species, 
      `Chlordanes BSAF (calc)` = Chlordane_calc,
      `Dieldrin BSAF (calc)` = Dieldrin_calc,
      `DDTs BSAF (calc)` = DDT_calc,
      `PCBs BSAF (calc)` = PCB_calc,
      `Chlordanes Conc (ng/g)` = Chlordane_conc,
      `Dieldrin Conc (ng/g)` = Dieldrin_conc,
      `DDTs Conc (ng/g)` = DDT_conc,
      `PCBs Conc (ng/g)` = PCB_conc
    )
  
  return(out)
  
}
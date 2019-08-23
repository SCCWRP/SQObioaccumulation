# install.packages('tidyverse')
library(tidyverse)

source('R/funcs.R')

# data inputs
# tischmthr <- read.csv('data/tischmthr.csv', stringsAsFactors = FALSE)
# finalsiteassess <- read.csv('data/finalsiteassess.csv', stringsAsFactors = FALSE)
# indic_lookup <- read.csv('data/indic_lookup.csv', stringsAsFactors = FALSE)
biota <- read.csv('data/biota.csv', stringsAsFactors = FALSE)
constants <- read.csv('data/constants.csv', stringsAsFactors = FALSE)
contam <- read.csv('data/contam.csv', stringsAsFactors = FALSE)
biota_preyprop <- read.csv('data/preyprop.csv', stringsAsFactors = FALSE, row.names = 1)

# calculated contaminant inputs
contamcalc <- cntcalc(contam, constants)

btch <- bioaccum_batch(biota, contamcalc, biota_preyprop, constants)

# append to reactive output
cbiota <- btch$CBIOTA
bsaf <- btch$BSAF
contamcalc <- contamcalc


  
  
  
  

# bsaf results
bsaf <- reactive({
  
  # input
  bsaf <-res$bsaf
  
  # bsaf
  nms <- colnames(bsaf) 
  out <- bsaf %>% 
    data.frame %>% 
    rownames_to_column('species') 
  names(out) <- c('species', nms)
  
  return(out)
  
})

# concentration in biota
cbiota <- reactive({
  
  # input
  cbiota <- res$cbiota
  
  # cbiota
  nms <- colnames(cbiota) 
  out <- cbiota %>% 
    data.frame %>% 
    rownames_to_column('species') 
  names(out) <- c('species', nms)
  
  return(out)
  
})

# combine bsaf and cbiota output
rescmb <- reactive({
  
  # input
  cntbsaf <- input$cntbsaf
  bsaf <- bsaf()
  cbiota <- cbiota()
  
  req(cntbsaf)
  
  # bsaf
  bsaf <- bsaf %>% 
    select(species, !!cntbsaf) %>% 
    rename(
      BSAF = !!cntbsaf
    ) %>% 
    mutate(species = factor(species, levels = species))
  
  # cbiota
  cbiota <- cbiota %>% 
    select(species, !!cntbsaf) %>% 
    rename(
      `Tissue conc. (ng/g wet)` = !!cntbsaf
    ) %>% 
    mutate(species = factor(species, levels = species))
  
  out <- bsaf %>% 
    left_join(cbiota, by = 'species') %>% 
    gather('var', 'val', -species)
  
  return(out)
  
})

# summary of contaminants across indicator guilds
indic_sum <- reactive({
  
  # input
  cbiota <- cbiota()
  contamcalc <- res$contamcalc
  
  req(!is.null(res$contamcalc))
  
  out <- indic_sum_fun(cbiota, contamcalc)
  
  return(out)
  
})


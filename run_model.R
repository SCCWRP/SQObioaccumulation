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
cbiota <- btch$cbiota
bsaf <- btch$bsaf

# summary table
indic_sum_fun(cbiota, contamcalc)

# plot of bsaf, cbiota by specific contaminant
plo_bsaf(bsaf, cbiota, 'alphaChlordane')

# tabular summary of bsaf, cbiota by specific contaminant
tab_bsaf(bsaf, cbiota, 'alphaChlordane')

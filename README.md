
# README

[![Travis build
status](https://travis-ci.org/SCCWRP/SQObioaccumulation.svg?branch=master)](https://travis-ci.org/SCCWRP/SQObioaccumulation)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/SCCWRP/SQObioaccumulation?branch=master&svg=true)](https://ci.appveyor.com/project/SCCWRP/SQObioaccumulation)

Install the package as follows:

``` r
install.packages('devtools')
library(devtools)
install_github('SCCWRP/SQObioaccumulation')
library(SQObioaccumulation)
```

Run the bioaccumulation model with defaults:

``` r
# data inputs
data(biota)
data(constants)
data(contam)
data(mcsparms)

# calculated contaminant inputs
contamcalc <- cntcalc(contam, constants)

# run model
res <- bioaccum_batch(biota, contamcalc, constants)

# assign output to separate objects
cbiota <- res$cbiota
bsaf <- res$bsaf
```

Creat a summary table:

``` r
# summary table
indic_sum <- indic_sum_fun(cbiota, contamcalc)
indic_sum
```

    ## # A tibble: 9 x 9
    ##   Guild `Chlordanes BSA~ `Dieldrin BSAF ~ `DDTs BSAF (cal~ `PCBs BSAF (cal~
    ##   <chr>            <dbl>            <dbl>            <dbl>            <dbl>
    ## 1 indi~             5.30            2.10             11.2              9.13
    ## 2 indi~             5.45            1.99             12.4             10.0 
    ## 3 indi~             5.26            2.25              9.55             7.93
    ## 4 indi~             6.33            3.22             13.0             11.1 
    ## 5 indi~             2.28            1.21              4.43             3.84
    ## 6 indi~             5.42            3.89              7.30             6.64
    ## 7 indi~             1.40            0.924             2.05             1.84
    ## 8 indi~             4.15            3.71              3.98             4.03
    ## 9 indi~             4.75            1.70             11.6              9.32
    ## # ... with 4 more variables: `Chlordanes Conc (ng/g)` <dbl>, `Dieldrin
    ## #   Conc (ng/g)` <dbl>, `DDTs Conc (ng/g)` <dbl>, `PCBs Conc (ng/g)` <dbl>

Plot BSAF and tissue concentration estimates for a selected contaminant:

``` r
# plot of bsaf, cbiota by specific contaminant
plo_bsaf(bsaf, cbiota, 'alphaChlordane')
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Make a table of BSAF and tissue concentration estimates for a selected
contaminant:

``` r
# tabular summary of bsaf, cbiota by specific contaminant
tab_bsaf(bsaf, cbiota, 'alphaChlordane')
```

    ## # A tibble: 2 x 28
    ##   Output Sediment Phytoplankton `Submerged Macr~ Zooplankton
    ##   <fct>     <dbl>         <dbl>            <dbl>       <dbl>
    ## 1 Tissu~      0.5         0.177            0.166       0.302
    ## 2 BSAF        1           0.355            0.331       0.605
    ## # ... with 23 more variables: `Small polychaete (e.g., Harmothoe
    ## #   imbricata)` <dbl>, `Large polychaete (e.g., Neanthes)` <dbl>,
    ## #   Amphipod <dbl>, Cumacean <dbl>, Mysid <dbl>, `Bivalve mollusk` <dbl>,
    ## #   `Decapod crab` <dbl>, `Crangon shrimp` <dbl>, `Forage fish -
    ## #   herbivore` <dbl>, `Forage fish - planktivore` <dbl>, `Forage fish -
    ## #   mixed diet i` <dbl>, `Forage fish - mixed diet ii` <dbl>, `Forage fish
    ## #   - primarily benthivore` <dbl>, `Forage fish - benthivore` <dbl>,
    ## #   indic1 <dbl>, indic2 <dbl>, indic3 <dbl>, indic4 <dbl>, indic5 <dbl>,
    ## #   indic6 <dbl>, indic7 <dbl>, indic8 <dbl>, indic9 <dbl>

Run Monte Carlo simulations (MCS) with results from bioaccumulation
model and additional inputs:

``` r
mcsres <- mcs_fun(1000, indic_sum, mcsparms, constants)
```

Summarize MCS results:

``` r
mcs_sum_fun(mcsres)
```

    ## # A tibble: 4 x 12
    ## # Groups:   Compound [4]
    ##   Compound    `0%`    `1%`   `5%`  `10%`  `25%` `50%` `75%` `90%` `95%`
    ##   <chr>      <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 Chlorda~ 0.147   0.208   0.287  0.347  0.501  0.755 1.14  1.58   1.93
    ## 2 DDT      0.151   0.210   0.385  0.524  0.909  1.59  2.74  4.70   6.54
    ## 3 Dieldrin 0.656   0.823   1.23   1.46   2.03   2.85  3.90  5.45   6.71
    ## 4 PCB      0.00417 0.00792 0.0223 0.0320 0.0636 0.142 0.322 0.760  1.15
    ## # ... with 2 more variables: `99%` <dbl>, `100%` <dbl>

Plot cumulative distribution curves for MCS:

``` r
mcs_plo(mcsres, xmax = 3)
```

<img src="README_files/figure-gfm/unnamed-chunk-9-1.png" width="80%" style="display: block; margin: auto;" />

Get overall SQO assessment:

``` r
wgtavg <- wgt_avg_fun(mcsparms)
sqo_sum_fun(wgtavg, mcsres, constants)
```

    ## # A tibble: 4 x 9
    ##   Compound `Observed tissu~ `Chemical expos~ `Estimated tiss~
    ##   <chr>               <dbl> <chr>                       <dbl>
    ## 1 Chlorda~             2.28 Very Low                    1.72 
    ## 2 DDT                  4.85 Very Low                    7.69 
    ## 3 Dieldrin             0.25 Very Low                    0.711
    ## 4 PCB                 36.5  Moderate                    5.17 
    ## # ... with 5 more variables: `Site linkage 25%` <dbl>, `Site linkage
    ## #   50%` <dbl>, `Site linkage 75%` <dbl>, `Site linkage category` <chr>,
    ## #   `Site assessment category` <chr>

### Metadata
Resources: <a href="https://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/1000_SQOHumanHealthFramework.pdf">report</a><br>
Contact: <a href="https://www.sccwrp.org/about/staff/ashley-parks/?pdf=4384">Ashley Parks</a>, <a href="stevenbay78@gmail.com">Steve Bay</a><br>

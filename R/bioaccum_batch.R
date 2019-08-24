#' Run bioaccumulation model function FoodWeb_SQO in batch for all species, contaminants
#'
#' @param biota input biological data
#' @param contam input contaminant data
#' @param biota_preyprop input prey proportions
#' @param constants input constants
#'
#' @import tibble
#' 
#' @importFrom magrittr "%>%"
bioaccum_batch <- function(biota, contam, biota_preyprop, constants){
  
  ### HAVE HARD-CODED THE PARAMETER NAMES
  theParamNames <- c("k1","k2", "GR", "Gv", "Gd", "Gf", "vlg", "vcg", "vng", "vwg", "kgb", 
                     "ke", "kd", "Ew", "Ed", "phi", "cpw", "assimEff_1", "assimEff_2", "assimEff_3");
  
  nspecies = nrow(biota_preyprop); 
  nprey = ncol(biota_preyprop);
  ncontam = length(contam$Chem); 
  nparams = length(theParamNames);
  
  #dietary composition (lipid, organic carbon, organic matter, water)
  biota_vld  = as.matrix(biota_preyprop) %*% biota$lipid;
  biota_vcd  = as.matrix(biota_preyprop) %*% biota$nloc;
  biota_vnd  = as.matrix(biota_preyprop) %*% biota$nlom;
  biota_wc = (1 - biota$lipid - biota$nloc - biota$nlom); #wc = (1 -lipid -nloc -nlom)
  #respiratory uptake of overlying water
  biota_mo = 1 - biota$mp;
  
  ## ASSIGN CONSTANT PARAMETERS from input data##
  EdA = constants$Value[which(constants$Constant=='EdA')];
  EdB = constants$Value[which(constants$Constant=='EdB')];
  xdoc = constants$Value[which(constants$Constant=='xdoc')];
  ddoc = constants$Value[which(constants$Constant=='ddoc')];
  xpoc = constants$Value[which(constants$Constant=='xpoc')];
  dpoc = constants$Value[which(constants$Constant=='dpoc')];
  alphapoc = constants$Value[which(constants$Constant=='alphapoc')];
  alphadoc = constants$Value[which(constants$Constant=='alphadoc')];
  ocsed = constants$Value[which(constants$Constant=='ocsed')]; 
  ds = constants$Value[which(constants$Constant=='ds')];
  A = constants$Value[which(constants$Constant=='A')];
  B = constants$Value[which(constants$Constant=='B')];
  T = constants$Value[which(constants$Constant=='T')];
  Cox = constants$Value[which(constants$Constant=='Cox')];
  vss = constants$Value[which(constants$Constant=='vss')];
  scav = constants$Value[which(constants$Constant=='scav')];
  
  #some descriptions of parameters
  ## that are used in model
  #cbiota = 'contaminant concentration in biota';
  #cprey  = 'contaminant concentration in prey diet';
  #k1 = 'aqueous uptake rate constant';
  #k2 = 'elimination rate constant';
  #GR = 'growth rate coefficient for fish or invertebrates; growth rate for phytoplankton';
  #Gv = 'gill ventilation rate';
  #Gd = 'feeding rate';
  #Gf = 'fecal egestion rate';
  #vlg = 'lipid fraction of gut';
  #vcg = 'nloc fraction of gut';
  #vng = 'nlom fraction of gut';
  #vwg = 'water fraction of gut';
  #kgb = 'gut-biota partition coefficient';
  #ke =  'fecal egestion rate constant (1/d)';
  #kd =  'dietary uptake rate constant';
  #Ew =  'gill chemical uptake efficiency';
  #Ed =  'dietary chemical transfer efficiency (also called gut uptake efficiency)';
  #phi = 'freely dissolved contaminant fraction in overlying water column';
  #cpw = 'contaminant concentration in porewater';
  #assimEff_1  = 'assimilation efficiency for lipid';
  #assimEff_2  = 'assimilation efficiency for nlom or nloc';
  #assimEff_3  = 'assimilation efficiency for water';
  
  ##Initialize output variables data objects (matrices and an array)
  #matrix that stores organism contaminant concentrations
  CBIOTA <- matrix(numeric(nspecies*ncontam), nrow=nspecies,dimnames=list(biota$Biota,contam$Chem));
  #contaminant concentrations of organisms' prey
  CPREY  <- matrix(numeric(nspecies*ncontam), nrow=nspecies,dimnames=list(biota$Biota,contam$Chem));
  #calculated biota sediment accumulation factor (BSAF)
  BSAF <-  matrix(numeric(nspecies*ncontam), nrow=nspecies,dimnames=list(biota$Biota,contam$Chem));
  #array stores each parameter specific to species and contaminant type
  PARAMS <- array(NA, c(nspecies,ncontam,nparams),dimnames=list(biota$Biota,contam$Chem,theParamNames));
  
  #PART 1
  ####Deterministic version###########################################################3
  
  #N.B. the use of for loops is rather inefficient computationally but was easy to code
  #and moderately easy to debug.  Perhaps there is a more efficient way.
  
  # loops over ALL contaminants
  # wrap the loop execution in withProgress
  for (icontam in 1:ncontam){
    
    # loop over species
    for (ispecies in 1:nspecies) {
      
      # # log
      # txt <- paste0('Contaminant', icontam, 'of', ncontam, '\n\tspecies', ispecies, 'of', nspecies)
      
      ## ASSIGN VARIABLES ##
      csed <- contam$cs_ng.g[icontam];  #sediment contaminant concentration 
      logkow_tempcor <- contam$logkow_tempcor[icontam]; #octanol-water partitioning coefficient (affects contaminant partitioning)
      log_KowTS <- contam$log_KowTS[icontam]; #temperature and salinity corrected octanol water partitioning
      Koc <- alphapoc*10^log_KowTS; # from Gobas and Arnot, 2010 Supp. pg. 5 (cit. Seth et at., 1999)
      
      ## Multiple water conc options ##
      cwater <- contam$free_cd_ng.ml[icontam]
      cpw <- contam$calc_cp_ng.ml[icontam]
      
      #other parameters, most reviewed above; described in Gobas and Arnot 2010
      beta <- contam$beta[icontam]; 
      betap <- contam$betap[icontam]; 
      kM <- contam$kM[icontam];
      phi <- contam$phi[icontam];
      taxa <- biota$taxa[ispecies];
      lipid <- biota$lipid[ispecies];
      nloc <- biota$nloc[ispecies];
      nlom <- biota$nlom[ispecies];
      wc <- biota_wc[ispecies]; 
      mp <- biota$mp[ispecies];
      mo <- biota_mo[ispecies];
      Wb <- biota$Wb[ispecies];
      GR <- biota$GR[ispecies];
      preyprop <- biota_preyprop[ispecies,]; # a vector of prey proportions
      cbiota <- CBIOTA[,icontam];
      vld <- biota_vld[ispecies]; # lipid proportions in diet
      vcd <- biota_preyprop[ispecies, 1] * ocsed + biota_vcd[ispecies]; # nonlipid O.C. proportions in diet
      vnd <- biota_vnd[ispecies]; # nonlipid O.M. proportions in diet
      vwd <- 1 - (vld + vcd + vnd); # water proportions in diet
      assimEff_1 <- biota$assimEff_1[ispecies]; # assimilation efficiency for lipid
      assimEff_2 <- biota$assimEff_2[ispecies]; # assimilation efficiency for non-lipid crabon (nloc) or org matter (nlom)
      assimEff_3 <- biota$assimEff_3[ispecies]; # assimilation efficiency for water
      
      ###########################
      ### CALL FOOD WEB MODEL ###
      
      # if(icontam == 1 & ispecies == 10) browser()
      
      #call the function, using all the various model input parameters
      Results <- FoodWeb_SQO(NumSim=1, csed, cwater, cpw, log_KowTS, logkow_tempcor, EdA, EdB, xdoc, ddoc, xpoc, dpoc, alphapoc, 
                             alphadoc, ocsed, ds, taxa, A, B, T, lipid, nloc, nlom, wc, beta, betap, mo, mp, phi, kM, Wb, Cox, 
                             vss, scav, preyprop, cbiota, vld, vcd, vnd, vwd, GR, assimEff_1, assimEff_2, assimEff_3) 
      
      #extract the biota and prey contaminant information
      CBIOTA[ispecies,icontam] <- Results$cbiota;
      CPREY[ispecies,icontam] <- Results$cprey;
      # check column headers on Results match dimnames on PARAMS 
      for (o in 1:(length(Results)-2)) 
      { if (dimnames(PARAMS)[[3]][o]==colnames(Results)[o+2]) {PARAMS[ispecies,icontam,o] <- Results[[o+2]]} else { print('Error in ordering of parameters')} 
      };  
    } #end of loop over species
  } #end of loop over contaminants
  
  
  ### CALCULATE BSAF (simple measure of biota versus sediment contamination) ###
  for (i in 1:nspecies) {
    BSAF[i,] <- CBIOTA[i,]/contam$cs_ng.g
  }
  
  ## minor formatting of bsaf and cbiota
  
  # bsaf
  nms <- colnames(BSAF) 
  bsaf <- BSAF %>% 
    data.frame %>% 
    rownames_to_column('species') 
  names(bsaf) <- c('species', nms)
  
  # cbiota
  nms <- colnames(CBIOTA) 
  cbiota <- CBIOTA %>% 
    data.frame %>% 
    rownames_to_column('species') 
  names(cbiota) <- c('species', nms)
  
  out <- list(bsaf = bsaf, cbiota = cbiota)
  
  return(out)
  
}

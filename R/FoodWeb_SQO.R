#' Core food web model
#'
#' @param NumSim required input
#' @param csed required input
#' @param cwater required input 
#' @param cpw required input
#' @param log_KowTS required input
#' @param logkow_tempcor required input
#' @param EdA required input
#' @param EdB required input
#' @param xdoc required input
#' @param ddoc required input
#' @param xpoc required input
#' @param dpoc required input
#' @param alphapoc required input
#' @param alphadoc required input
#' @param ocsed required input
#' @param ds required input
#' @param taxa required input
#' @param A required input
#' @param B required input
#' @param T required input
#' @param lipid required input
#' @param nloc required input
#' @param nlom required input
#' @param wc required input
#' @param beta required input
#' @param betap required input
#' @param mo required input
#' @param mp required input
#' @param phi required input
#' @param kM required input
#' @param Wb required input
#' @param Cox required input
#' @param vss required input
#' @param scav required input
#' @param preyprop required input
#' @param cbiota required input
#' @param vld required input
#' @param vcd required input
#' @param vnd required input
#' @param vwd required input
#' @param GR required input
#' @param assimEff_1 required input
#' @param assimEff_2 required input
#' @param assimEff_3 required input
#'
#' @export
#'
FoodWeb_SQO <- function(NumSim, csed, cwater, cpw, log_KowTS, logkow_tempcor, EdA, EdB, xdoc, 
                        ddoc, xpoc, dpoc, alphapoc, alphadoc, ocsed, ds, taxa, A, B, T, lipid, nloc, 
                        nlom, wc, beta, betap, mo, mp, phi, kM, Wb, Cox, vss, scav, preyprop, cbiota, 
                        vld, vcd, vnd, vwd, GR, assimEff_1, assimEff_2, assimEff_3){ #updated 6.30.2010
  
  ## initialize the output variables
  temp <- data.frame(cbiota = NA, cprey = NA, k1=NA, k2=NA, GR=NA, Gv=NA, Gd=NA, 
                     Gf=NA, vlg=NA, vcg=NA, vng=NA, vwg=NA, kgb=NA, ke=NA, kd=NA, Ew=NA, Ed=NA, phi=NA, 
                     cpw=NA, assimEff_1=NA, assimEff_2=NA, assimEff_3=NA);
  Output <- temp[rep(1:nrow(temp), times=NumSim), ]
  
  lden <- 0.9; # lipid density - added 6.30.2010
  
  Kow <- 10^logkow_tempcor # kow temp corrected, arithmetic
  KowTS <- 10^log_KowTS # Kow temp and salinity corrected, arithmetic
  
  ##
  ## Calculate contaminant dependant parameters
  ##
  #Ew <- 1/(1.85+1.55/KowTS); # gill chemical uptake efficiency
  Ew <- 1/(1.85+155/KowTS); # gill chemical uptake efficiency #fixed 6.30.2010
  Ed <- 1/(EdA*Kow + EdB); # dietary chemical transfer efficiency (gut uptake efficiency)
  
  ## assign contaminant specific outputs
  Output$Ew <- Ew;
  Output$Ed <- Ed;
  Output$phi <- phi;
  Output$cpw <- cpw;
  
  ##
  ## Now proceed with the biota-specific calculations.
  ##
  ## Calculations are taxa dependant. The main IF loop here represents that dependancy.
  
  if (taxa==0) {Output$cbiota <- csed}; # This is a dummy taxa used by sediment
  
  if (taxa==1)
  {
    # plants (i.e., phytoplankton) only have selected parameters
    
    # aqueous uptake rate constant different for phytoplankton
    k1 <- ((A + (B/KowTS))^-1);
    
    # biota-water partition coefficient for phytoplankton; i.e., the bioconcentration factor
    #kbw <- (lipid*KowTS + nlom*betap*KowTS + wc); 		
    kbw <- (lipid*KowTS/lden + nloc*betap*KowTS + nlom*beta*KowTS + wc); #updated 6.30.2010
    
    # elimination rate constant for phytoplankton ******AM changed from kbwp********
    k2 <- k1/kbw; 
    kG <- GR
    # calculate concentration in organism
    # based on uptake/loss pseudoequilibrium
    # main model equation (only includes plant parameters)
    
    Output$cbiota <- (k1 * cwater) / (k2 + kG)
    # assign output variables
    Output$k1 <- k1;
    Output$k2 <- k2;
    #		Output$kG <- kG;
  }
  
  if (taxa>=2)  # ML- CHANGED THIS FROM taxa>1 b/c code was trying calculated fecal egestion rate for taxa 1.1 which is sea lettuce (ulva)!!! However, should check that this doesn't cause other problems. ###
  {
    if (taxa==2|taxa==3|taxa==4|taxa==5|taxa==6) 
    {
      Gv <- (1400*Wb^0.65)/Cox;  # Calculate Gv - gill ventilation rate
      k1 <- Ew*Gv/Wb;  # k1 = aqueous uptake rate constant (L/kg*d)            
      #kbw <- (lipid*KowTS + nlom*beta*KowTS + wc);  # kbw = biota-water partition coefficient
      kbw <- (lipid*KowTS/lden + nloc*betap*KowTS + nlom*beta*KowTS + wc);  #updated 6/30/2010 & fixed 7/1/2010
      
      k2 <- k1/kbw;	 # k2 = gill elimination rate constant
      
      # calculate feeding rate
      if (taxa==3|taxa==5) 
      {
        
        Gd <- 0.022 * (Wb^0.85) * exp(0.06*T); # non filter feeders
        
      }
      if(taxa==6){
        
        Gd <- (0.022 * (Wb^0.85) * exp(0.06*T) + Gv*vss*scav) / 2; # forage fish - herbivore
        
      }
      if (taxa==2|taxa==4) 
      {
        Gd <- Gv*vss*scav; # filter feeders (kg/d)
      }
      
      kd <- Ed*Gd/Wb;	# calculate kd dietary uptake rate constant (kg/kg*d)
    }
    
    kG <- GR * Wb^-0.2;
    #-----------------
    
    # calculate fecal egestion rate (kg/d)
    Gf <- Gd*((1-assimEff_1)*vld + (1-assimEff_2)*(vnd+vcd) + (1-assimEff_3)*vwd); #updated 6/30/2010
    
    # lipid (vlg), nloc (vcg), nlom (vng), and water (vwg) fraction of gut contents (kg/kg), respectively
    ## Updated 6/30/2010 ##
    vlg <- (1-assimEff_1)*vld/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd);
    vng <- (1-assimEff_2)*vnd/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd); 
    vcg <- (1-assimEff_2)*vcd/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd);
    vwg <- (1-assimEff_3)*vwd/((1-assimEff_1)*vld+(1-assimEff_2)*(vnd+vcd)+(1-assimEff_3)*vwd);
    
    kgb <- (vlg*Kow/lden + vcg*betap*Kow + vng*beta*Kow + vwg)/
      (lipid*Kow/lden + nlom*beta*Kow + wc); # gut-biota partition coefficient (unitless) - updated 6/30/2010
    
    ke <- Gf*Ed*kgb/Wb; # fecal egestion rate constant (1/d)
    
    # calculate total concentration available from prey for dietary uptake
    Output$cprey <- as.vector(as.matrix(preyprop) %*% cbiota);
    
    # finally, calculate concentration in organism based on uptake/loss pseudo-equilibrium.
    # Output$cbiota <- (k1*(mo*phi*cwater + mp*cpw)+kd*Output$cprey[1])/(k2 + ke + kG + kM);
    
    Output$cbiota <- (k1*((1 - mp)*cwater + mp*cpw)+kd*Output$cprey[1])/(k2 + ke + kG + kM);	
    
    # assign results to out variable
    Output$k1  <- k1;
    Output$k2  <- k2;
    #Output$kG  <- kG;
    Output$GR  <- GR;
    Output$Gv  <- Gv;
    Output$Gd  <- Gd;
    Output$Gf  <- Gf;
    Output$assimEff_1 <- assimEff_1;
    Output$assimEff_2 <- assimEff_2;	
    Output$assimEff_3 <- assimEff_3;	
    Output$vlg <- vlg;
    Output$vcg <- vcg;
    Output$vng <- vng;
    Output$vwg <- vwg;
    Output$kgb <- kgb;
    Output$ke  <- ke;
    Output$kd  <- kd;
    
  }
  return(Output)
  
}


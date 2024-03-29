#
# Authors: Jake Diamond and Enrico Bertuzzo
# Purpose: A model for dissolved oxygen, CO2, pH, and carbonate sytems
# Date: 2023 February 17
# 
source(file.path("src", "functions_for_model.R"))
# Overall model function ----------------------------------------------------------
model <- function(time, states, pars,
                  light_forcing = FALSE, temp_forcing = FALSE){
  with(as.list(c(states, pars)), {
    ################## Forcing functions if used ##########################
    if (temp_forcing) {
      temp <- temp_signal(time)
    } else { 
      temp <- temp
    }
    if (light_forcing) {
      # Photosynthetically reactive radiation; PAR (umol/m2/s)
      par  <- light_signal(time)
      # mean daily PAR
      meanpar <- mean(light_signal(0:23))
    } else {
      doy <- doy + floor(time/24)
      # mean PAR at the time step (umol/m2/s)
      par <- par_fun(time, day = doy, latitude = latitude)
      # mean daily PAR
      meanpar <- mean(par_fun(0:23, day = doy, latitude = latitude))
    }

    ################ Molecular weights ################
    MW_CO2    <- 44.01 # (g/mol) molecular weight CO2
    MW_O2     <- 31.998 # (g/mol) molecular weight O2
    MW_CaCO3  <- 100.0869 # (g/mol) molecular weight CaCO3
    MW_Ca     <- 40.078 # (g/mol) molecular weight CaCO3
    
    ################ Calculate internal parameters ################
    I <- 1.6E-5 * cond # Ionic strength (mol/m3)
    S <- 54 * I # Salinity (ppt)
    
    # Water density (kg/m3)
    rho_w <- rhow(temp)

    # Calculate Schmidt numbers for CO2 and O2 (Rik Wanninkhof L&O methods 2014)
    Sc_CO2 <- Sc("CO2", temp)
    Sc_O2  <- Sc("O2", temp)
    
    # Calculate gas exchange coefficients for CO2 and O2 (1/d)
    K_CO2 <- K600 * (600 / Sc_CO2)^0.5
    K_O2  <- K600 * (600 / Sc_O2)^0.5
    
    # Calculate gas exchange velocity (m/h)
    k_CO2 <- K_CO2 * z / (24) 
    k_O2  <- K_O2 * z / (24) 
    
    # Temperature in Kelvin
    TK <- temp + 273.15
    
    # -log(monovalent act. coef)
    gamma <- act(TK, cond) 
    
    # Calculate Henry's constant K0 for CO2 (mol/kg/atm); Weiss 1974 Mar. Chem.
    # using freshwater version from Plummer and Busenberg (1982)
    K_henry <- KH(TK)
    
    # Calculate saturation value for O2
    O2_sat <- O2_sat(temp, p_atm) / MW_O2 # (mol/m3)
    
    # Calculate saturation value for CO2
    CO2_sat <- CO2_atm * 1e-6 * p_atm * K_henry * rho_w # (mol/m3)
    
    #################### Calculate carbonate equilibrium #####################
    # compute pH, HCO3, and CO2 from alkalinity and DIC
    carb_eq <- carb(TK = TK, AT = ALK, pH = 8, cond = cond, TC = DIC)
    
    # Get the equilibrium values of the carbonate system
    pH    <- carb_eq$pH
    CO2   <- carb_eq$CO2
    HCO3  <- carb_eq$HCO3
    CO3   <- carb_eq$CO3
    Ca    <- HCO3 / 2 #assuming 1:2 relationship in calcium carbonate system

    # Chemical enhancement of fCO2
    if ( alf ) {
    alpha_enh <- chem_enh(temp = temp, pH = pH, KCO2 = K_CO2, d = z)
    } else { alpha_enh <- 1 }
    
    # Calcite saturation index (unitless), concentrations in mol/L = log(IAP/Ksp)
    SI <- log10(((10^(4 * -gamma) * Ca / 1000)) * (10^(4 * -gamma) * CO3 / 1000 ) / 
                  Ksp(TK, sal = S))
    
    ############################ Calculate fluxes ###########################
    # First calculate GPP as linear function of mean PAR (mol O2/m2/hr)
    gpp <- (gpp_mean / 24 / MW_O2) * (par / meanpar)
    
    # Ecosystem respiration (mol O2/m2/hr)
    er <- (er_mean / 24) / MW_O2
    
    # Gas exchange (mol/m2/del_t) from river perspective (+ is into river)
    fO2  <- k_O2 * (O2_sat - O2) * open
    fCO2 <- k_CO2 * (CO2_sat - CO2) * alpha_enh * open
    
    # exchange with lateral inflow  (mol/m2/hr)
    LO2  <- qL_w * (O2_l - O2)
    LDIC  <- qL_w * (DIC_l - DIC)
    LALK <- qL_w * (ALK_l - ALK)

    # Calcite dissolution (mol CaCO3/m3/hr); account for stoichiometry in
    # CaCO3 + CO2 + H2O <--> Ca2++ + 2HCO3- 
    # The change in DIC is one mol per mol CaCO3 precipitated/dissolved
    # This, calc, returns the rate of CaCO3 dissolution
    calc <- calc_rate(temp, I, SI, pH, CO2 / 1000, CALC, S)
    # calc <- 0
    ################### Difference equations ################################
    # Rates of change (mol/m3/hr)
    dO2   <- ((gpp - er + fO2 + LO2) / z)
    dDIC  <- ((-gpp * (1 / PQ) + er * RQ + fCO2 + LDIC) / z + calc)
    dCALC <- -calc
    dALK  <- 2 * calc + LALK / z # dALK is two mol per mol CaCO3 ppt/dissolved
    
    ################### Non-CO2 supported GPP ################################
    # If the change in DIC is greater than available [CO2] when GPP is using DIC
    # Then autotrophs are using HCO3
    if ( -dDIC > CO2 && gpp > 0 ) {
      dHCO3 <- CO2 + dDIC
      # How much of GPP is supported by HCO3?
      GPP_HCO3 <- -dHCO3 / (gpp / z) * 100
    } else {
      GPP_HCO3 <- 0 # GPP supported by HCO3, initial value
    }
    
    ################### Model output ################################
    # Differences
    diffs <- c(dO2 = dO2, dDIC = dDIC, dALK = dALK, dCALC = dCALC)
    
    # Fluxes
    fluxes <- c(GPP  = gpp,
                ER   = -er,
                NEP  = (gpp - er),
                fO2  = -fO2,
                fCO2 = -fCO2,
                LO2  = LO2,
                LDIC = LDIC,
                LALK = LALK,
                calc = -calc)

    # Other variables of interest
    CO2sat <- CO2 / CO2_sat * 100
    O2sat  <- O2 / O2_sat * 100
    exCO2  <- CO2 - CO2_sat
    exO2   <- O2 - O2_sat
    # exDIC  <- DIC - DIC_sat

    # Output of model
    return(list(diffs, fluxes, PAR = par, pH = pH, 
                CO2 = CO2, HCO3 = HCO3, CO3 = CO3,
                CO2sat = CO2sat, O2sat = O2sat, temp = temp, GPP_HCO3 = GPP_HCO3,
                exCO2 = exCO2, exO2 = exO2, enh = alpha_enh, K_H = K_henry,
                SI = SI)) #exDIC = exDIC, 
    return(list(diffs))
  })
}  # end of model

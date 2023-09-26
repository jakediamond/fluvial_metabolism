#
# Authors: Jake Diamond and Enrico Bertuzzo
# Purpose: To model dissolved oxygen, CO2, pH, and carbonate sytems
# Date: 2023 February 17
# 

# Load libraries
library(lubridate)
library(tidyverse)
library(deSolve)

# Get other functions -----------------------------------------------------
source(file.path("src", "LSI_function.R"))
source(file.path("src", "co2_chem_enh_function.R"))
source(file.path("src", "results_helper_functions.R"))
source(file.path("src", "calcite_kinetics.R"))
source(file.path("src", "functions.R"))

# Get simulation times ----------------------------------------------------
# Get simulation times
del_t    <- 0.5               # time step (h)
days     <- 2                 # number of days to simulate
sim_time <- 24 * days / del_t # simulation time (-)
times    <- seq(0, sim_time)  # sequence of in-model times to simulate

# Define model parameters --------------------------------------------------------
parms <- c(
  # geographic parameters
  latitude = 47.8,
  doy = 190, # day of the year
  
  # Hydraulic parameters
  z = 1, # (m) average depth
  
  # Physicochemical parameters
  CO2_atm = 400,   # (uatm) atmospheric pressure of CO2
  p_atm   = 1,     # (atm) barometric atmospheric pressure
  temp    = 20,    # (deg C) temperature of water
  pH      = 7,     # pH of the water
  Alk     = 1.5,   # (mol m-3)total measured alkalinity 
  Ca      = 0.75,  # (mol m-3) calcium concentration
  cond    = 300,   # (uS/cm) specific conductance
  
  # Reaction parameters
  gpp_mean = 10,  # (g O2/m2/d) daily mean GPP
  K600     = 1,   # (1/d) gas exchange coefficient at Schmidt number 600
  er_mean  = 10,  # (g O2/m2/d) daily mean ER
  PQ       = 1,   # photosynthetic quotient relating mols CO2 to produce mols O2
  RQ       = 1    # respiratory quotient relating mols O2 to produce mols CO2
)


# Determine model forcing if desired --------------------------------------
# Stream temperature
temp_signal <- approxfun(x = data$time_hr / del_t, 
                         y = data$temp, 
                         method = "linear", rule = 2)
# Stream light
light_signal <- approxfun(x = data$time_hr / del_t, 
                          y = data$light, 
                          method = "linear", rule = 2)
# Overall model function ---------------------------------------------------
model <- function(time, states, parms, 
                  light_forcing = FALSE, temp_forcing = FALSE){
  with(as.list(parms), {
    
    # Unpack states
    O2 <- states[1]
    DIC <- states[2]
    ALK <- states[3]
    
    ################## Forcing functions if used ##########################
    if (temp_forcing) {
      temp <- temp_signal(time)
    } else { 
      temp <- temp
    }
    if (light_forcing) {
      # Photosynthetically reactive radiation; PAR (umol/m2/s)
      par  <- light_signal(time)
      # Total sum of PAR for that day
      sumpar <- sum(light_signal(0:((24/del_t) - 1)) * del_t)
    } else {
      # mean PAR at the time step (umol/m2/s)
      par <- par_fun(time * del_t, day = doy, latitude = latitude)
      # Total sum of PAR for that day
      sumpar <- sum(par_fun(0:((24/del_t) - 1), 
                            day = doy, latitude = latitude) * del_t)
    }
    
    ################ Molecular weights ################
    MW_CO2    <- 44.01 # (g/mol) molecular weight CO2
    MW_O2     <- 31.998 # (g/mol) molecular weight O2
    MW_CaCO3  <- 100.0869 # (g/mol) molecular weight CaCO3
    MW_Ca     <- 40.078 # (g/mol) molecular weight CaCO3
    
    ################ Calculate internal parameters ################
    # Water density (kg/m3)
    rho_w <- rhow(temp)
    
    # Useful conversions for later
    ALK_molkg <- ALK / rho_w # Total alkalinity from mol/m3 to mol/kg
    DIC_molkg <- DIC / rho_w # DIC in mol/kg
    
    # Calculate Schmidt numbers for CO2 and O2 (Rik Wanninkhof L&O methods 2014)
    Sc_CO2 <- Sc("CO2", temp)
    Sc_O2  <- Sc("O2", temp)
    
    # Calculate gas exchange coefficients for CO2 and O2 (1/d)
    K_CO2 <- K600 / ((600 / Sc_CO2)^-0.5)
    K_O2  <- K600 / ((600 / Sc_O2)^-0.5)
    
    # Calculate gas exchange velocity (m/del_t)
    k_CO2 <- K_CO2 * z / (24 / del_t) 
    k_O2  <- K_O2 * z / (24 / del_t) 
    
    # Temperature in Kelvin
    TK <- temp + 273.15
    
    # Calculate Henry's constant K0 for CO2 (mol/kg/atm); Weiss 1974 Mar. Chem.
    # using version from Plummer and Busenberg (1982)
    K_henry <- henry(TK)
    
    # Calculate saturation value for O2
    O2_sat <- O2_sat(temp, p_atm) / MW_O2 # (mol/m3)
    
    # Calculate saturation value for CO2
    CO2_sat <- CO2_atm * 1e-6 * p_atm * K_henry * rho_w # (mol/m3)
    
    # Calculate DIC equilibrium with atmosphere
    DIC_sat <- seacarb::carb(flag = 24, var1 = CO2_atm, var2 = ALK_molkg,
                             S = 0, T = temp, warn = "n")$DIC * rho_w # (mol/m3)
    
    #################### Calculate carbonate equilibrium #####################
    # compute pH, HCO3, and CO2 from alkalinity and DIC using seacarb
    carb_eq <- seacarb::carb(flag = 15, var1 = ALK_molkg, var2 = DIC_molkg,
                             S = 0, T = temp, warn = "n")
    # Use initial conditions for pH and HCO3 if possible
    pH    <- carb_eq$pH
    CO2   <- carb_eq$CO2 * rho_w # (mol/m3)
    HCO3  <- carb_eq$HCO3 * rho_w # (mol/m3)
    
    # Chemical enhancement of fCO2
    alpha_enh <- chem_enh_fun(temp = temp, pH = pH, KCO2 = K_CO2, d = z)
    
    # Langelier Saturation Index (unitless)
    LSI <- LSI(temp = temp, pH = pH, HCO3 = HCO3 / 1000, cond = cond, Ca = Ca) 

    ############################ Calculate fluxes ###########################
    # First calculate GPP as linear function of total PAR (mol O2/m2/del_t)
    gpp <- (gpp_mean / MW_O2) * (par / sumpar)

    # Ecosystem respiration amount (mol O2/m2/del_t)
    er <- (er_mean / 24) / MW_O2
    
    # Gas exchange (mol/m2/del_t) from river perspective (+ is into river)
    fO2  <- k_O2 * (O2_sat - O2)
    fCO2 <- k_CO2 * (CO2_sat - CO2) * alpha_enh
    
    # Calcite precipitation (mol C/m2/del_t); account for stoichiometry in
    # Ca2++ + 2HCO3- <--> CaCO3 + CO2 + H2O
    # This is the moles of HCO3 lost from solution
    calc <- calc_rate_fun(temp, cond, pH, Ca, CO2 / 1000, HCO3 / 1000)
    calc_dic <- -calc
    
    ################### Difference equations ################################
    # Rates of change (mol/m3)
    dO2   <- ((gpp - er + fO2) / z) * del_t 
    dDIC  <- ((-gpp * PQ + er * RQ + fCO2 - calc_dic) / z) * del_t
    
    # Determine if Alkalinity changes 
    # If the change in DIC is greater than available [CO2]
    # Then autotrophs are using HCO3, so ALK is reduced
    if (-dDIC > CO2) {
      dALK <- CO2 + dDIC
    } else {
      dALK <- -calc_dic
    } 
    
    # Differences
    diffs <- c(dO2 = dO2, dDIC = dDIC, dALK = dALK)
    
    # Fluxes
    fluxes <- c(GPP  = gpp  * del_t, 
                ER   = er   * del_t, 
                fO2  = fO2  * del_t,
                fCO2 = fCO2 * del_t,
                calc = calc_dic * del_t)
    
    # Other variables of interest
    CO2sat <- CO2 / CO2_sat * 100
    O2sat  <- O2 / O2_sat * 100
    exCO2  <- CO2 - CO2_sat
    exO2   <- O2 - O2_sat
    exDIC  <- DIC - DIC_sat
    
    # Output of model
    return(list(diffs, fluxes, PAR = par, pH = pH, CO2 = CO2, HCO3 = HCO3,
                CO2sat = CO2sat, O2sat = O2sat, temp = temp,
                exCO2 = exCO2, exO2 = exO2,
                exDIC = exDIC, LSI = LSI))
  })
}  # end of model

# Initial conditions ------------------------------------------------------
# Equilibrium based on pH and ALK from parameters
carb_ini <- seacarb::carb(flag = 8,
                          var1 = parms["pH"],
                          var2 = parms["Alk"] / 1000,
                          T = parms["temp"],
                          S = 0,
                          warn = "n")

# Stream DO concentration and carbonate system
O2_ini  <- 9 / 32               # mol/m3, ~ atmospheric eq
DIC_ini <- carb_ini$DIC * 1000  # mol/m3  
ALK_ini <- parms["Alk"]         # mol/m3

# Initial conditions vector
yini <- c(O2 = O2_ini,
          DIC = DIC_ini,
          ALK = ALK_ini)

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
ode_out <- ode(y = yini,
               times = times,
               func = model,
               parms = parms,
               light_forcing = F,
               temp_forcing = F)

# Examine the model output ------------------------------------------------
# Get into good formats
output <- ode_to_df(ode_out)

# Check out some plots
plot_fun(output, c("exCO2", "exO2"), "biplot")
plot_fun(output, c("exDIC", "exO2"), "biplot")
plot_fun(output, c("exDIC", "exO2"), "timeseries")
plot_fun(output, "calc")
plot_fun(output, c("ALK", "DIC"))
plot_fun(output, "pH")
plot_fun(output, c("fCO2", "calc"))
plot_fun(output, "ALK")
plot_fun(output, "CO2")
plot_fun(output, "exCO2")
plot_fun(output, "LSI")
plot_fun(output, "CO2sat")
plot_fun(output, "temp")
plot_fun(output, "O2")
plot_fun(output, "PAR")



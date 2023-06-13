#
# Authors: Jake Diamond and Enrico Bertuzzo
# Purpose: Run the metabolism/carbonate model
# Date: 2023 February 17
# 

# Load libraries
library(tidyverse)
library(deSolve)

# Get model functions -----------------------------------------------------
source(file.path("src", "functions_for_model.R"))
source(file.path("src", "results_helper_functions.R"))
source(file.path("src", "carbonate_model.R"))

# Get simulation times ----------------------------------------------------
del_t    <- 1                # time step (h)
days     <- 10               # number of days to simulate
sim_time <- (24 * days ) - 1 # simulation time (h)
times    <- seq(0, sim_time, # sequence of in-model times to simulate
                by = del_t)  

# Define model parameters --------------------------------------------------------
pars <- c(
  # geographic parameters
  latitude = 47.8,
  doy = 236, # day of the year
  
  # Hydraulic parameters
  Q    = 33.8, #(m3 s-1) average discharge
  z    = 1,    # (m) average depth
  w    = 200, # (m) average width of stream
  QL_f = 0.2,  # (-) fraction of discharge from lateral inflow 
  qL   = 3E-4, # (m2/s) average lateral inflow per unit length
  
  # Physicochemical parameters
  CO2_atm = 411,   # (uatm) atmospheric pressure of CO2
  p_atm   = 1.012, # (atm) barometric atmospheric pressure
  temp    = 20,    # (deg C) temperature of water
  pH      = 6.75,  # pH of the water
  Alk     = 1.5,     # (mol m-3) total measured alkalinity 
  Ca      = 1.4,   # (mol m-3) calcium concentration
  CaCO3   = 0,     # (mol m-3) calcite concentration
  cond    = 317,   # (uS/cm) specific conductance
  ALK_l   = 1.5,     # (mol m-3) ALK concentration of lateral inflow
  O2_l    = 0.1,   # (mol m-3) O2 concentration of lateral inflow
  DIC_l   = 2,     # (mol m-3) DIC concentration of lateral inflow
  
  # Reaction parameters
  gpp_mean = 12,  # (g O2/m2/d) daily mean GPP
  K600     = 3,  # (1/d) gas exchange coefficient at Schmidt number 600
  er_mean  = 10,  # (g O2/m2/d) daily mean ER
  PQ       = 1,  # photosynthetic quotient relating mols O2 to produce mols CO2
  RQ       = 1,  # respiratory quotient relating mols O2 to produce mols CO2
  alf      = 1,  # whether (1) or not (0) to include chemical fco2 enhance
  open     = 1,  # whether system is open (1) or closed (0)
  qL_known = 1   # whether lateral discharge is known (1) or not (0)
)

# Determine model forcing if desired --------------------------------------
# mods_out <- readRDS(file.path("results", "model_comparison_results_florida.RDS"))
# x <- pluck(mods_out, 3, 2)
# a <- pull(x[1:192, "time_hr"] / 0.25)
# a <- a[seq(1, length(a), 4)]
# a <- 0:47
# b <- pull(x[1:192, "temp"])
# b <- b[seq(1, length(b), 4)]
# c <- pull(x[1:192, "light"])
# c <- c[seq(1, length(c), 4)]
# 
# # Stream temperature
# temp_signal <- approxfun(x = a,
#                          y = b,
#                          method = "linear", f = 0.5)
# # Stream light
# light_signal <- approxfun(x = a,
#                           y = c,
#                           method = "linear", f = 0.5)


# Initial conditions ------------------------------------------------------
# Equilibrium based on pH and ALK from parameters
# carb_ini <- carb(TK = pars["temp"] + 273.15, AT = pars["Alk"], 
#                  pH = pars["pH"], cond = pars["cond"])

# Equilibrium with atmosphere and alkalinity
carb_ini <- seacarb::carb(flag = 24,
                          var1 = pars["CO2_atm"], var2 = pars["Alk"] / rhow(pars["temp"]),
                          pHscale = "F", k1k2 = "m06",
                          S = pars["cond"] * 1.6E-5 * 54,
                          T = pars["temp"], warn = "n")


# Initial O2 concentration and carbonate system (mol/m3)
O2_ini   <- as.numeric(O2_sat(pars["temp"])) / 31.998  
DIC_ini  <- as.numeric(carb_ini$DIC * rhow(pars["temp"]))
ALK_ini  <- as.numeric(carb_ini$ALK * rhow(pars["temp"]))
CALC_ini <- as.numeric(pars["CaCO3"]) 

# Initial conditions vector
yini <- c(O2 = O2_ini,
          DIC = DIC_ini,
          ALK = ALK_ini,
          CALC = CALC_ini)

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
# ode_out <- ode(y = yini,
#                times = times,
#                func = model,
#                parms = pars,
#                # light_forcing = TRUE, temp_forcing = TRUE
# )
# # 
# # # # Get into good formats
# output <- ode_to_df(ode_out)
# # 
# # # # Check out some plots
# plot_fun(output, "DIC")

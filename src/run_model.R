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
del_t    <- 0.25                # time step (h)
days     <- 2               # number of days to simulate
sim_time <- (24 * days ) - 1 # simulation time (h)
times    <- seq(0, sim_time, # sequence of in-model times to simulate
                by = del_t)  

# Define model parameters --------------------------------------------------------
pars <- c(
  # geographic parameters
  latitude = 46.59,
  doy = 266, # day of the year
  
  # Hydraulic parameters
  Q    = 11.95, # (m3 s-1) average discharge
  z    = 0.53,    # (m) average depth
  w    = 29,  # (m) average width of stream
  qL   = 3.5E-5, # (m2/s) average lateral inflow per unit length
  
  # Physicochemical parameters
  CO2_atm = 411,   # (uatm) atmospheric pressure of CO2
  p_atm   = 1.012, # (atm) barometric atmospheric pressure
  temp    = 20,    # (deg C) temperature of water
  pH      = 7.24,  # pH of the water
  Alk     = 1.03,   # (mol m-3) total measured alkalinity 
  Ca      = 1.4,   # (mol m-3) calcium concentration
  CaCO3   = 0.2,     # (mol m-3) calcite concentration
  cond    = 430,   # (uS/cm) specific conductance
  ALK_l   = 3.5,   # (mol m-3) ALK concentration of lateral inflow
  O2_l    = 0.2,   # (mol m-3) O2 concentration of lateral inflow
  DIC_l   = 3,     # (mol m-3) DIC concentration of lateral inflow
  
  # Reaction parameters
  gpp_mean = 12, # (g O2/m2/d) daily mean GPP
  K600     = 2.6,  # (1/d) gas exchange coefficient at Schmidt number 600
  er_mean  = 12, # (g O2/m2/d) daily mean ER
  PQ       = 1.2,  # photosynthetic quotient relating mols O2 to produce mols CO2
  RQ       = 1.5,  # respiratory quotient relating mols O2 to produce mols CO2
  alf      = 1,  # whether (1) or not (0) to include chemical fco2 enhance
  open     = 1  # whether system is open (1) or closed (0)
)

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
O2_ini   <- 0.29#as.numeric(O2_sat(pars["temp"])) / 31.998  
DIC_ini  <- 1.16#as.numeric(carb_ini$DIC * rhow(pars["temp"]))
ALK_ini  <- 1.03#as.numeric(carb_ini$ALK * rhow(pars["temp"]))
CALC_ini <- as.numeric(pars["CaCO3"]) 

# Initial conditions vector
yini <- c(O2 = O2_ini,
          DIC = DIC_ini,
          ALK = ALK_ini,
          CALC = CALC_ini)

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
ode_out <- ode(y = yini,
               times = times,
               func = model,
               parms = pars,
               light_forcing = FALSE, temp_forcing = FALSE
)

# Get into good formats
output <- ode_to_df(ode_out)

# # # Check out some plots
# plot_fun(output)
# plot_fun(output, "pH")
# plot_fun(output, "DIC")
# # plot_fun(output, "HCO3")
# plot_fun(output, "GPP_HCO3")
# # plot_fun(output, "ALK")
# # plot_fun(output, "SI")
# plot_fun(output, "NEP")
# # plot_fun(output, c("DIC", "exO2"), "biplot")
# # plot_fun(output, c("HCO3", "exO2"), "biplot")
# # 
# # 
# # tail(output, 24) |>
# #   # filter(exO2 > 0) |>
# #   # summarize(mean(GPP_HCO3))
# #   ggplot(aes(x = HCO3, y = exO2)) +
# #   geom_point()
# 
# 
# summary(lm(exO2~DIC, data = tail(output, 24)))
# summary(lm(exO2~HCO3, data = tail(output, 24)))
# summary(lm(exO2~HCO3, data = tail(output, 24) |>
#              filter(exO2 < 0 )))
# 

ggplot() +
  geom_line(data = output,
            aes(x = time,
            y = O2)) +
  geom_point(data = rename(x, time = time_hr),
             aes(x = time,
                 y = O2))

ggplot() +
  geom_line(data = output,
            aes(x = time,
                y = CO2)) +
  geom_point(data = rename(x, time = time_hr),
             aes(x = time,
                 y = CO2))


plot(temp_signal(0:28))
plot(light_signal(0:48))
plot(output$PAR)


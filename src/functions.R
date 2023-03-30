# -------------------------------------
# Author: Jake Diamond
# Purpose: support functions for parameter estimates in metabolism model
# Date: 2023-03-01
# -------------------------------------

# Density of water function
rhow <- function (temp = 25) {
  # Martin and McCutcheon (1999)
  999.842594 + 6.793952*10^(-2)*temp - 9.095290*10^(-3)*temp^2 + 
    1.001685*10^(-4)*temp^3 - 1.120083*10^(-6)*temp^4 + 6.536335e-9*temp^5
}

# Schmidt number function
Sc <- function (gas = "CO2", temp = 25) {
  co <- switch(gas,
    "O2" =  c(1923.6, -125.06, 4.3773, -0.085681, 0.00070284),
    "CO2" = c(1745.1, -125.34, 4.8055, -0.101150, 0.00086842)
  )
  co[1] + co[2] * temp + co[3] * temp^2 + co[4] * temp^3 + co[5] * temp^4
}


# Function for Henry's constant at given temperature and pressure (mol/kg/atm)
henry <- function (TK = 298, press = 1) {
  # TK = temperature in Kelvin, 
  # press = pressure in atm
  # using version from Plummer and Busenberg (1982)
  K0 <- 10^(108.3865 + 0.01985076*TK - 6919.53/TK - 
              40.45154*log10(TK) + 669365/TK^2)
  # Correct for pressure; Weiss 1974 Marine Chemistry
  R <- 82.05736 # Gas constant, [cm3 * atm / (K * mol)]
  vbarCO2 <- 32.3 #partial molal volume of CO2 in solution (cm^3/mol)
  KH <- K0 * exp(((1 - press) * vbarCO2)/(R * TK))
  return(KH)
}

# O2 saturation function
O2_sat <- function (temp = 25, press = 1) {
  # Conversions
  mlL_mgL   <- 1.42905 # O2 mL/L to mg/L, per USGS memo 2011.03
  mmHg_atm  <- 760 # mmHg to atm
  
  # Calculate saturated concentrations of O2 (mg/L) using Garcia-Benson
  # vapor pressure of water (mmHg)
  u <- 10^(8.10765 - 1750.286 / (235 + temp)) 
  #press correction (=1 at 1atm) USGS memos 81.11 and 81.15
  press_corr <- (press * mmHg_atm - u) / (760 - u) 
  Ts <- log((298.15 - temp)/(273.15 + temp)) # scaled temperature
  lnC <- 2.00907 + 3.22014 * Ts + 4.0501 * Ts^2 + 4.94457 * 
    Ts^3 + -0.256847 * Ts^4 + 3.88767 * Ts^5
  o2.sat <- exp(lnC) #O2 saturation (mL/L)
  o2_sat <- o2.sat * mlL_mgL * press_corr # (mg/L)
  return(o2_sat)
}

# Photosynthetically active radiation function (umol/m2/s)
par_fun <- function (hour, day = 180, latitude = 0, max.insolation = 2326) {
  declin <- (23.45 * sin((360 / 365) * (284 + day))) * pi / 180
  hour.angle <- (360/24) * (hour - 12) * pi / 180
  lat <- latitude * pi / 180
  zenith <- acos(sin(lat) * sin(declin) + 
                   cos(lat) * cos(declin) * cos(hour.angle))
  insolation <- max.insolation * cos(zenith)
  insolation <- pmax(insolation, 0)
}

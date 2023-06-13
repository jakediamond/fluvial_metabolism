# -------------------------------------
# Author: Jake Diamond
# Purpose: Functions to solve the carbonate system for use in model
# Date: 2023-04-10
# -------------------------------------

# Helper functions --------------------------------------------------------
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
  declin <- (23.45 * sin((2*pi / 365) * (284 + day))) * pi / 180
  hour.angle <- (360/24) * (hour - 12) * pi / 180
  lat <- latitude * pi / 180
  zenith <- acos(sin(lat) * sin(declin) + 
                   cos(lat) * cos(declin) * cos(hour.angle))
  insolation <- max.insolation * cos(zenith)
  insolation <- pmax(insolation, 0)
}

# Density of water function kg/m3
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

# -log10(activity coeff.), gamma, function
# TK is temperature in Kelvin
# cond is specific conductance in uS/cm
# I is ionic strength in mol/L, estimated from cond unless known
act <- function(TK, cond, I = NULL){
  if(is.null(I)) {
    I <- 1.6*10^-5 * cond
  }
  E <- 60954 / (TK + 116) - 68.937 #dielectric constant
  # Davies approximation for activity coefficient
  gamma <- 1.82 * 10^6 * (E * TK)^-1.5 * ((sqrt(I) / (1 + sqrt(I))) - 0.3 * I)
  return(gamma)
}

# Calculate first freshwater thermodynamic equilibrium constant
# [H+][HCO3-]/[CO2]
# at infinite dilution via Plummer and Busenberg 1982
K1 <- function(TK) {
  pK1t <- 356.3094 + 0.06091964*TK - 21834.37 / TK - 126.8339 * 
    log10(TK) + 1684915/TK^2
  K1t <- 10^-pK1t
  return(K1t)
}

# Calculate second freshwater thermodynamic equilibrium constant
# [[H+][CO3--]/[HCO3-]
# at infinite dilution via Plummer and Busenberg 1982
K2 <- function(TK){
  pK2t <- 107.8871 + 0.03252849*TK - 5151.79 / TK - 38.92561 * 
    log10(TK) + 563713.9/TK^2
  K2t <- 10^-pK2t
  return(K2t)
}

# Kw function, thermodynamic dissociation constant for water
Kw <- function(TK) {
  # Kw is dissociation constant for water
  pKw <- 4471 / TK + 0.01706 * TK -6.0875
  Kw <- 10^-pKw
}

# Function for Henry's constant at given temperature and pressure (mol/kg/atm)
KH <- function (TK = 298, press = 1) {
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

# Apparent, or stoichiometric, solubility constant for calcite
# Mucci 1983
# Ca2++ + CO3-- = CaCO3 (s)
Ksp <- function (TK = 298.15, sal = 0) {
  # this is the -log(thermodynamic solubility constant) as a function of temp. 
  # in distilled water according to Plummer and Busenberg 1982
  pKsp_0 <- -171.9065 - 0.077993 * TK + 2839.319/TK + 71.595 * log10(TK)
  # These are how salinity affects the constant according to 
  # Mucci 1983, Table 7
  B <- +(-0.77712 + 0.0028426 * TK + 178.34/TK) * sqrt(sal)
  C <- -0.07711 * sal + 0.0041249 * sal^1.5
  log10Kspc <- pKsp_0 + B + C
  Kspc <- 10^(log10Kspc)
}

# Carbonate system function -----------------------------------------------
# Calculate carbonate system based on alkalinity and pH
# TK = temperature in Kelvin
# AT = total alkalinity (mol/m3)
# pH = pH on free scale
# cond = specific conductance at 25degC (uS/cm)
carb <- function(TK, AT, pH, cond) {
  # Calculate all the parameters
  AT <- AT / rhow(TK-273.15) # mol/kg
  gamma <- act(TK, cond) # -log(monovalent act. coef)
  Kh <- KH(TK) #henry's constant uncorrected for salinity
  I <- 1.6*10^-5 * cond # ionic strength, cond in uS/cm
  S <- 53.974*I #salinity from ionic strength, estimated
  KHa <- Kh + (0.023517 - 0.023656 * TK/100 + 0.0047036 * TK/100 * TK/100)*S #apparent Henry constant, Weiss 1974
  aH <- 10 ^ -gamma # activity coefficient for H+
  aOH <- 10 ^ -gamma # activity coefficient for OH-
  aHCO3 <- 10 ^ -gamma # activity coefficient for HCO3-
  aCO3 <- 10 ^ (4 * -gamma) # activity coefficient for CO32-
  H <- 10^(-pH) #hydrogen ion conc.
  K1t <- K1(TK) #thermodynamic diss. coeff.
  K2t <- K2(TK) #thermodynamic diss. coeff.
  KWa <- Kw(TK) / (aH * aOH) # apparent dissociation coefficient 
  K1a <- K1(TK) / (aH * aHCO3) # apparent dissociation coefficient 
  K2a <- K2(TK) / (aH * aCO3 / aHCO3) # apparent dissociation coefficient 
  OH <- KWa / H
  
  # Solve the carbonate system, from Stumm and Morgan 
  alpha1 <- (H * K1a)/(H^2 + K1a * H + K1a * K2a)  # HCO3 ionization fraction 
  alpha2 <- (K1a * K2a) / (H^2 + K1a * H + K1a * K2a) # CO3 ionization fraction
  CT  <- (AT- OH + H) / (alpha1 + 2*alpha2) #total carbon, DIC
  CO2 <- CT * (H ^ 2) / (H ^ 2 + K1a * H + K1a * K2a)
  HCO3 <- CO2 * K1a / H
  CO3 <- HCO3 * K2a / H
  
  # Get concentrations in mol/m3
  CO2 <- CO2 * rhow(TK - 273.15)
  HCO3 <- HCO3 * rhow(TK - 273.15)
  CO3 <- CO3 * rhow(TK - 273.15)
  DIC <- CT * rhow(TK - 273.15)
  
  #CO2 in mol/kg
  CO2_molkg <- (CT * (H^2) * 10^(6*gamma)) / ((H^2 * 10^(6*gamma)) + 
                                           (K1t * H * 10^(4*gamma)) + (K1t * K2t))
  
  # #Uses Henry's Law constant and converts from atm to uatm (KH in fugacity (mol-atm / kg-
  # #soln))
  # pCO2 <- (CO2 / KH) * 1000000
  
  #Uses the apparent Henry's Law constant and converts from atm to uatm (KHa in fugacity
  #(mol-atm / kg-soln))
  pCO2 <- (CO2_molkg / KHa) * 1000000
  # Return the carbonate system
  data.frame(CO2 = CO2, pCO2 = pCO2, DIC = DIC, HCO3 = HCO3, CO3 = CO3)
}

# Chemical enhancement of fCO2 --------------------------------------------
# temperature dependent rate constants for CO2 hydration/dehydration
rate_const <- function(TK, h, sal = 0){
  # The combined rate constant, from Hoover and Berkshire 1969 is:
  # r = kCO2 + kOH * [OH-]
  # and [OH-] = Kw / [H+]
  # following Johnson 1982 L&O
  # CO2 + H2O = H2CO3
  lnkCO2 <- 1246.98 + 0 * sal - 6.19 * 10^4 / TK - 183 * log(TK)
  kCO2 <- exp(lnkCO2) #1/s
  # CO2 + OH- = HCO3-
  # This empirical estimate has the dissociation constant for water built in (Kw)
  lnkOH.Kw <- -930.13 + 0.11 * sal + 3.1 * 10^4 / TK + 140.9 * log(TK)
  kOH.Kw <- exp(lnkOH.Kw) #mol/L/s
  # Combined rate constant
  r <- kCO2 + kOH.Kw / h
  return(r)
}

# temperature dependent CO2 diffusion coefficient (m2/s)
D <- function(TK){
  # following Zeebe 2011 Geochemica et Cosmochimica Acta
  D <- 0.0000000146836 * ((TK / 217.2056) - 1)^1.997
  return(D)
}

# chemical enhancement function based on: 
# temperature (degC)
# pH
# CO2 gas exchange coefficient, KCO2 (1/d)
# depth, d (m)
chem_enh <- function(temp, pH, KCO2, d, cond = 300) {
  TK <- temp + 273.15 #temperature in Kelvin
  pfm <- act(TK, cond) # -log10(monovalent activity)
  h <- 10^(pfm - pH) # activity of hydrogen ion
  k_cms <- KCO2 * d * 100 / 86400 # CO2 gas exchange constant in cm/s
  # following Hoover and Berkshire 1969 and Wanninkhoff and Knox 1996
  # simplified variable from Hoover and Berkshire 1969
  tau <- 1 + h^2 / (K1(TK) * K2(TK) + K1(TK) * h)
  r <- rate_const(TK, h) #1/s
  D <- D(TK) * 1e4 # diffusivity of CO2 (cm2/s)
  Q <- sqrt(r * tau / D) # 1/cm
  Z <- D / k_cms # cm, boundary layer thickness
  alpha <- tau / ((tau - 1) + tanh(Q * Z) / (Q * Z)) # dimensionless
  return(alpha)
}

# Calcite kinetics --------------------------------------------------------
calc_rate <- function (temp, I, SI, pH, CO2, CaCO3) {
  # This equation is from Plummer et al. 1978 10.2475/ajs.278.2.179
  # Also known as the "PWP" equation
  # R = k1 * aH + k2 * aH2CO3star + k3 * aH20 - k4 * aCa * aHCO3
  # where ki are kinetic coefficients that are a function of temperature and ax 
  # is the activity of ion x. This results in a a surface-area specific rate of mols/cm2
  # Assume aH2CO3star is equal to activity of CO2
  # Assume activity of CO2 and H2O are 1
  # temperature in Kelvin
  TK <- temp + 273.15
  
  # # Calculate LSI
  # LSI <- LSI(temp, pH, cond, Ca, HCO3)
  # ionic strength (mol/m3), cond in uS/cm
  # I <- 1.6*10^-5 * cond 
  # activity coefficient of uncharged species (e.g., CO2, H2CO3)
  a_uncharge <- 10^(0.1 * I)
  
  # calculate kinetic rate constants
  # Equations 5-9 in Plummer et al. 1978
  logk1 <- 0.198 - 444 / TK 
  logk2 <- 2.84 - 2177 / TK
  if (temp < 25) {
    logk3 <- -5.86 - 317 / TK
  } else {
    logk3 <- -1.10 - 1737 / TK
  }
  k1 <- 10^logk1 # cm/s
  k2 <- 10^logk2 # cm/s
  k3 <- 10^logk3 # mmol/cm2/s
  # calcite solubility constant
  Ksp <- Ksp(TK)
  # K2 carbonate equilibrium constant
  K2 <- K2(TK)
  # -log10(monovalent activity coefficient)
  pfm <- act(TK, I=I)
  fm <- 10^-pfm
  # hydrogen ion activity
  aH <- 10^(pfm - pH)
  # carbonic acid activity
  aH2CO3star <- CO2 * a_uncharge
  # water activity = 1
  # aH2O <- 1
  # calcium activity
  # aCa <- Ca * fm^2
  # HCO3 activity
  # aHCO3 <- HCO3 * fm
  # Equation 25
  k4 <- (K2 / Ksp) * (k1 + (1 / aH) * (k2 * aH2CO3star + k3))
  # Forward reaction rate (dissolution) (mmol/cm2/s)
  rf <- k1 * aH + k2 * aH2CO3star + k3
  # Use simplifying assumptions; Parkhurst and Appelo 1999 (PHREEQc V2 p. 43)
  rk <- rf * (1 - 10^(2/3 * SI))
  rk <- rk / 1000 * 3600 #(mol/cm2/hr), per area of calcite
  
  # The overall rate for a kinetic reaction of minerals and other solids is
  # Rk = rk * (A0/V) * (m / m0)^n
  # rk is the specific rate (mol/m2/s)
  # A0 is the initial surface area of the solid (m2)
  # V is the amount of solution (kg)
  # mo is the initial moles of solid
  # m is the moles of the solid at a given time
  # (m/mo)^n is a factor to account for changes in A0/V during
  # dissolution and aging of the solid
  # For uniformly dissolving spheres and cubes n = 2/3
  
  # Specific surface area of calcite (cm^2/mol), this should be a parameter
  # as it has a strong control on the rate, but is set for now
  sa <- 5000 #PWP 1978
  
  # Initial amount of CaCO3 (mol/m3)
  # CaCO3_0 <- 0
  # Calculate kinetics only if there is some CaCO3 to dissolve,
  # or if LSI allows precipitation
  if (CaCO3 <= 0 & SI < 0){
    rate <- 0
  } else {
    # GETTING RID OF THIS FOR NOW, BUT PERHAPS SHOULD BE INCLUDED
    # # get the actual area (cm2) based on the mol/m3 of calcite 
    # if (CaCO3_0 > 0) {
    #   area <- sa * CaCO3_0 * (CaCO3 / CaCO3_0)^(2/3)
    # } else {
    area <- sa * CaCO3
    # }
    
    # Get rate in mol/m3 by accounting for the surface area of calcite
    rate <- area * rk # (mol/m3/hr)
  }
}






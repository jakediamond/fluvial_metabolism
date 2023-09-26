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
carbonate <- function(TK, AT, DIC, cond = 300, pH=8) {
  
  AT <- AT / rhow(TK - 273.15) # mol/kg
  DIC <- DIC / rhow(TK - 273.15) # mol/kg
  gamma <- act(TK, cond) # -log(monovalent act. coef)
  Kh <- KH(TK) #henry's constant uncorrected for salinity
  I <- 1.6 * 10^-5 * cond # ionic strength, cond in uS/cm
  S <- 53.974 * I #salinity from ionic strength, estimated
  aH <- 10 ^ -gamma # activity coefficient for H+
  aOH <- 10 ^ -gamma # activity coefficient for OH-
  aHCO3 <- 10 ^ -gamma # activity coefficient for HCO3-
  aCO3 <- 10 ^ (4 * -gamma) # activity coefficient for CO32-
  KWa <- Kw(TK) / (aH * aOH) # apparent dissociation coefficient
  K1a <- K1(TK) / (aH * aHCO3) # apparent dissociation coefficient
  K2a <- K2(TK) / (aH * aCO3 / aHCO3) # apparent dissociation coefficient
  KHa <- Kh + (0.023517 - 0.023656 * TK / 100 + 0.0047036 * TK / 100 * TK / 100) * S #apparent Henry constant, Weiss 1974
  
  # INTERNAL VARIABLES USED TO CALCULATE pCO2 and pH
  #       H      #  current concentration of hydrogen ion (mol/kg)
  #       CA     #  carbonate alkalinity (eq / l)
  #       CO2aq  #  concentration of aqueous CO2
  #       diff.H  #  difference in successive estimates of H+ (mol/kg)
  #       tiny.diff.H   (1e-15) test for convergence 
  #       a      #  first term in quadratic for eq 12
  #       b      #  second term in quadratic for eq 12
  #       c      # third term in quadratic for eq 12

  # Iterate for H and CA by repeated solution of eqs 13 and 12
  H <- 10.^(-pH)                 # initial guess from arg list      
  diff.H <- H     
  tiny.diff.H <- 1.e-15 
  
  iter <- 0
  
  while (diff.H > tiny.diff.H) {     # iterate until H converges
    
    H.old <- H                      # remember old value of H
    
    # solve Tans' equation 13 for carbonate alkalinity from TA
    CA <- AT  
    
    # solve quadratic for H (Tans' equation 12)
    a <- CA
    b <- K1a * (CA - DIC)
    c <- K1a * K2a * (CA - 2 * DIC)
    H <- (-b + sqrt(b^2 - 4. * a * c) ) / (2. * a)  
    
    # How different is new estimate from previous one?
    diff.H <- abs(H - H.old)
    iter <- iter + 1
    
  }
  
  # Solve the carbonate system, from Stumm and Morgan
  alpha1 <- (H * K1a) / (H^2 + K1a * H + K1a * K2a)  # HCO3 ionization fraction
  alpha2 <- (K1a * K2a) / (H^2 + K1a * H + K1a * K2a) # CO3 ionization fraction
  CO2 <- DIC * (H ^ 2) / (H ^ 2 + K1a * H + K1a * K2a)
  HCO3 <- CO2 * K1a / H
  CO3 <- HCO3 * K2a / H
  
  # Get concentrations in mol/m3
  CO2_mM <- CO2 * rhow(TK - 273.15)
  pH <- -log10(H)
  HCO3 <- HCO3 * rhow(TK - 273.15)
  CO3 <- CO3 * rhow(TK - 273.15)
  DIC <- DIC * rhow(TK - 273.15)
  AT <- AT * rhow(TK - 273.15)
  
  # Uses the apparent Henry's Law constant and converts from atm to uatm
  pCO2 <- (CO2 / KHa) * 1000000
  
  # Return the carbonate system
  data.frame(CO2 = CO2_mM, pCO2 = pCO2, pH = pH, DIC = DIC, HCO3 = HCO3, CO3 = CO3, AT = AT)
}

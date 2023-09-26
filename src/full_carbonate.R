carb2 <- function(TK, AT, pH, cond, TC = NULL) {
  
  # Calculate all the parameters
  AT <- AT / rhow(TK - 273.15) # mol/kg
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
  
  # Calculate pH if DIC is given
  if(!is.null(TC)) {
    TC <- TC / rhow(TK - 273.15) # mol/kg
    # Iterate for H and CA by repeated solution
    H <- 10 ^ (-pH)  # initial guess from arg list      
    delH <- H     
    tol <- 1.e-15 
    
    iter <- 0
    
    while (delH > tol) {     # iterate until H converges
      
      H_old <- H                      # remember old value of H
      
      # solve for carbonate alkalinity from TA
      CA <- AT  
      
      # solve quadratic for H
      a <- CA
      b <- K1a * (CA - TC)
      c <- K1a * K2a * (CA - 2 * TC)
      H <- (-b + sqrt(b^2 - 4 * a * c) ) / (2 * a)  
      
      # How different is new estimate from previous one?
      delH <- abs(H - H_old)
      iter <- iter + 1
      
    }
    pH <- -log10(H)
    CT <- TC
  } else {
    H <- 10^(-pH) #hydrogen ion conc.
  }
  
  OH <- KWa / H # hydroxide ion conc.
  
  # Solve the carbonate system, from Stumm and Morgan, (mol/kg)
  alpha1 <- (H * K1a) / (H^2 + K1a * H + K1a * K2a)  # HCO3 ionization fraction
  alpha2 <- (K1a * K2a) / (H^2 + K1a * H + K1a * K2a) # CO3 ionization fraction
  if(is.null(TC)) {
    CT  <- (AT - OH + H) / (alpha1 + 2 * alpha2) #total carbon, DIC
  }
  CO2 <- CT * (H ^ 2) / (H ^ 2 + K1a * H + K1a * K2a)
  HCO3 <- CO2 * K1a / H
  CO3 <- HCO3 * K2a / H
  
  # Get concentrations in mol/m3
  CO2_mM <- CO2 * rhow(TK - 273.15)
  HCO3 <- HCO3 * rhow(TK - 273.15)
  CO3 <- CO3 * rhow(TK - 273.15)
  DIC <- CT * rhow(TK - 273.15)
  
  # Uses the apparent Henry's Law constant and converts from atm to uatm
  pCO2 <- (CO2 / KHa) * 1000000
  
  # Return the carbonate system
  data.frame(pH = pH, CO2 = CO2_mM, pCO2 = pCO2, DIC = DIC, HCO3 = HCO3, CO3 = CO3)
}


seacarb::carb(flag = 15, var1 = ALK_molkg, var2 = DIC_molkg,
              pHscale = "F",
              S = S, T = temp, warn = "n")

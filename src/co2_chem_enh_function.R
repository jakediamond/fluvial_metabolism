# temperature dependent dissociation constants for carbonic acid
diss_const_fun <- function(TK, h, sal = 0){
  K1 <- K1_fun(TK)
  K2 <- K2_fun(TK)
  # simplified variable from Hoover and Berkshire 1969
  tau <- 1 + h^2 / (K1 * K2 + K1 * h)
}

# temperature dependent rate constants for CO2 hydration/dehydration
rate_const_fun <- function(TK, h, sal = 0){
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
}

# temperature dependent CO2 diffusion coefficient (m2/s)
D_fun <- function(TK){
  # following Zeebe 2011 Geochemica et Cosmochimica Acta
  D <- 0.0000000146836 * ((TK / 217.2056) - 1)^1.997
}

# chemical enhancement function based on: 
# temperature (degC)
# pH
# CO2 gas exchange coefficient, KCO2 (1/d)
# depth, d (m)
chem_enh_fun <- function(temp, pH, KCO2, d, cond = 300) {
  TK <- temp + 273.15 #temperature in Kelvin
  pfm <- pfm_fun(TK, cond = 300) # -log10(monovalent activity)
  h <- 10^(pfm - pH) # activity of hydrogen ion
  k_cms <- KCO2 * d * 100 / 86400 # CO2 gas exchange constant in cm/s
  # following Hoover and Berkshire 1969 and Wanninkhoff and Knox 1996
  tau <- diss_const_fun(TK, h) # dimensionsless
  r <- rate_const_fun(TK, h) #1/s
  D <- D_fun(TK) * 1e4 # diffusivity of CO2 (cm2/s)
  Q <- sqrt(r * tau / D) # 1/cm
  Z <- D / k_cms # cm, boundary layer thickness
  alpha <- tau / ((tau - 1) + tanh(Q * Z) / (Q * Z)) # dimensionless
  return(alpha)
}

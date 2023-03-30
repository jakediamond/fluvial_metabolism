# -------------------------------------
# Author: Jake Diamond
# Purpose: Calculate the calcite saturation index 
# Date: 2023 February 2
# -------------------------------------

# Mucci 1983
# Ca2++ + CO3-- = CaCO3 (s)
# Apparent, or stoichiometric, solubility constant for calcite
Ksp_fun <- function (TK = 298.15, sal = 0) {
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

# K1 function, first dissociation constant
K1_fun <- function(TK, sal = 0){
  # following Millero et al. 2002 Deep-Sea Research
  # [H+][HCO3-]/[CO2]
  # pK1 = -8.712-9.46*0.001*sal+8.56*0.00001*sal^2+1355.1/TK+1.7976*log(TK)
  # Different calculation for freshater: American Public Health Association 2005  
  # via JAWWA 1990, 82(7) via Plummer and Busenberg 1982
  pK1 <- 356.3094 + 0.06091964*TK - 21834.37 / TK - 126.8339 * 
    log10(TK) + 1684915/TK^2
  K1 <- 10^-pK1
}

# K2 function, second dissociation constant
K2_fun <- function(TK, sal = 0){
  # [H+][CO3--]/[HCO3-]
  # pK2 = 17.0001-0.01259*sal-7.9334*0.00001*sal^2+936.291/TK-1.87354*log(TK)-
  #   2.61471*sal/TK+0.07479*sal^2/TK
  # Different calculation for freshwater
  # via Plummer and Busenberg 1982
  pK2 <- 107.8871 + 0.03252849*TK - 5151.79 / TK - 38.92561 * 
    log10(TK) + 563713.9/TK^2
  K2 <- 10^-pK2
}

# Kw function, dissociation constant for water
Kw_fun <- function(TK) {
  # Kw is dissociation constant for water
  pKw <- 4471 / TK + 0.01706 * TK -6.0875
  Kw <- 10^-pKw
}

# pfm function, fm is activity coefficient for monovalent species
pfm_fun <- function(TK, cond){
  I <- 1.6*10^-5 * cond # ionic strength, cond in uS/cm
  E <- 60954/(TK+116) - 68.937 #dielectric constant
  # Davies approximation for activity coefficient
  pfm <- 1.82 * 10^6 * (E*TK)^-1.5 * ((sqrt(I)/(1+sqrt(I))) - 0.3 * I)
}

# Calcite saturation function. This predicts whether water has a tendency to 
# precipitate or dissolved CaCO3. This version is the Langelier Saturation Index,
# so if it is greater than 0, precipitation is likely, and if less than 0
# dissolution is likely
LSI <- function (temp = 25, #degC
                       pH = 7, 
                       cond = 250, #uS/cm
                       Ca = 0.00075, #mol/L
                       HCO3 = 0.0015, #mol/L
                       sal = 0) #ppt
  {
  # American Public Health Association 2005
  # via JAWWA 1990, 82(7)
  # SI_calcite = pH - pHs; where pHs is the pH of the water when in equilibrium 
  # with calcium carbonate at the existing [Ca2++] and [HCO3-]
  # pHs = pK2 - pKs + p[Ca2++] + p[HCO3-] + 5pf_m
  # Where "p" refers to the -log10()
  
  # Some quick unit changes
  pCa    <- -log10(Ca)
  pHCO3  <- -log10(HCO3)
  TK     <- temp + 273.15 # temperature in kelvin
  
  # K2 is second dissociation constant for carbonic acid
  pK2 <- -log10(K2_fun(TK))
  
  # Ks is solubility product constant for CaCO3
  pKs <- -log10(Ksp_fun(TK, sal))
  
  # Kw is dissociation constant for water
  pKw <- -log10(Kw_fun(TK))
  
  # fm is activity coefficient for monovalent species
  pfm <- pfm_fun(TK, cond)
  
  # final calculations
  pHs <- pK2 - pKs + pCa + pHCO3 + 5*pfm
  SI <- pH - pHs
  return(SI)
}



# function for calcium carbonate precipitation potential (CCPP) [mg/L]
# This determines the actual amount of calcium carbonate that precipitates
# given ideal nucleation conditions and assuming equilibrium is achieved
# and can be used then to determine how much
# CO2 is released from this precipitation reaction
# Ca2++ + 2HCO3- <--> CaCO3 + CO2 + H2O
# from JAWWA 1990, 82(7) via Plummer and Busenberg 1982
# Only consideres HCO3-, CO3--, OH-, and H+ in calculating alkalinity, and 
# does not consider ion pairs. Therefore it tends to overestimate CaCO3 that can
# be precipitated and underestimates the amount that can be dissolved
CCPP_fun <- function(temp, pH, cond, alk, Ca){
  # First determine properties of water
  TK = temp + 273.15 # temperature in kelvin
  
  # monovalent activity coefficient
  pfm = pfm_fun(TK, cond)
  
  # Initial alkalinity (eq/L) from mmol/L
  Alk_i = alk / 1000
  
  # Initial calcium concentration
  MW_Ca = 40.078 #g/mol
  Ca_i = Ca / 1000 / MW_Ca #mol/L
  
  # initial hydrogen ion concentration (mol/L)
  H_i = 10^(pfm - pH)
  
  # get pK for equilibrium
  pK1 = pK1_fun(TK)
  pK2 = pK2_fun(TK)
  pKw = pKw_fun(TK)
  pKs = -log10(Ksp_fun(TK))
  
  # get apparent equilibrum constants
  K1 = 10^(2*pfm - pK1)
  K2 = 10^(4*pfm - pK2)
  Kw = 10^(2*pfm - pKw)
  Ks = 10^(8*pfm - pKs)
  
  # determine values for parameters, p, s, and t for initial conditions
  p = (2 * H_i + K1) / K1
  s = H_i - (Kw / H_i)
  t = (2 * K2 + H_i)/ H_i
  
  # calculate initial acidity
  Acy_i = ((Alk_i + s) / t) * p + s
  
  # Use iterative calculations to determine pH when equilibrium is established
  # Assume a pH for the equilibrated condition (pHeq) and calculate the [H+]
  # Start with pHeq = 7
  pH_start = 7
  eq_fun = function(pH_eq) {
    # Get equilbrium [H+]
    H_eq = 10^(pfm - pH_eq)
    # Get equilibrium parameters
    p_eq = (2 * H_eq + K1) / K1
    r_eq = (H_eq + 2 * K2) / K2
    s_eq = H_eq - (Kw / H_eq)
    t_eq = (2 * K2 + H_eq)/ H_eq
    # get RHS and LHS of eqn from Plummer and Busenberg 1982
    RHS = 2 * Ks * r_eq * p_eq / (t_eq * (Acy_i - s_eq)) - (t_eq * (Acy_i - s_eq)/
                                                              p_eq) + s_eq
    LHS = 2 * Ca_i - Alk_i
    # Minimization function, minimize sum of squared errors
    sum(RHS - LHS)^2
  }
  # optimize to find the equilibrium pH
  pH_eq = optim(par = pH_start, fn = eq_fun)$par
  # Get equilbrium [H+]
  H_eq = 10^(pfm - pH_eq)
  # Get equilibrium parameters
  p_eq = (2 * H_eq + K1) / K1
  r_eq = (H_eq + 2 * K2) / K2
  s_eq = H_eq - (Kw / H_eq)
  t_eq = (2 * K2 + H_eq)/ H_eq
  # This difference must be positive
  # if (Acy_i - s_eq < 0)
  #   warning("not a valid solution")
  
  # Determine equilibrium alkalinity
  Alk_eq = (t_eq / p_eq) * (Acy_i - s_eq) - s_eq
  
  # Determine CCPP (mg/L)
  CCPP = 50000 * (Alk_i - Alk_eq)
}

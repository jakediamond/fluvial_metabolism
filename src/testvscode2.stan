functions {
    // // estimate photosynthetically active radiation function (umol/m2/s)
    // real par_fun(real hour, real day, real latitude) {
    // real declin = (23.45 * sin((2 * pi() / 365) * (284 + day))) * pi() / 180;
    // real hour_angle = (360/24) * (hour - 12) * pi() / 180;
    // real lat = latitude * pi() / 180;
    // real zenith = acos(sin(lat) * sin(declin) + cos(lat) * cos(declin) * cos(hour_angle));
    // real insol = 2326 * cos(zenith);
    // real insolation = max(insol, 0);
    // return insolation;
    // }
    
    // calculate oxygen saturation
    real o2_sat(real temp) {
    // Conversions
    real ml_mg = 1.42905; // O2 mL/L to mg/L, per USGS memo 2011.03
    real mm_atm = 760;  // mmHg to atm
    real press_cor = 1; 
    // Calculate saturated concentrations of O2 (mg/L) using Garcia-Benson
    // vapor pressure of water (mmHg)
    real u = 10^(8.10765 - 1750.286 / (235 + temp));
    real ts = log((298.15 - temp) / (273.15 + temp)); //scaled temperature
    real ln_c = 2.00907 + 3.22014 * ts + 4.0501 * ts^2 + 4.94457 *
      ts^3 + -0.256847 * ts^4 + 3.88767 * ts^5;
    real eo2_sat = exp(ln_c); //O2 saturation (mL/L)
    real o2_sat = eo2_sat * ml_mg * press_cor; // (mg/L)
    return o2_sat;
    }
    
    // -log10(activity coefficient)
    real act(real TK, real cond) {
    real I = 1.6 * 10^-5 * cond;
    real E = 60954 / (TK + 116) - 68.937; //dielectric constant
    // Davies approximation for activity coefficient
    real gamma = 1.82 * 10^6 * (E * TK)^-1.5 * ((sqrt(I) / (1 + sqrt(I))) - 0.3 * I);
    return gamma;
    }
    
    // Calculate first freshwater thermodynamic equilibrium constant
    // [H+][HCO3-]/[CO2]
    // at infinite dilution via Plummer and Busenberg 1982
    real K1(real TK) {
    real pK1 = 356.3094 + 0.06091964 * TK - 21834.37 / TK - 126.8339 *
    log10(TK) + 1684915 / TK^2;
    real k1 = 10^-pK1;
    return k1;
    }
    
    // Calculate second freshwater thermodynamic equilibrium constant
    // [[H+][CO3--]/[HCO3-]
    // at infinite dilution via Plummer and Busenberg 1982
    real K2(real TK){
    real pK2 = 107.8871 + 0.03252849 * TK - 5151.79 / TK - 38.92561 *
    log10(TK) + 563713.9 / TK^2;
    real k2 = 10^-pK2;
    return k2;
    }
    
    // Kw function, thermodynamic dissociation constant for water
    real Kw(real TK) {
    // Kw is dissociation constant for water
    real pKw = 4471 / TK + 0.01706 * TK - 6.0875;
    real kw = 10^-pKw;
    return kw;
    }
    
    // Function for Henry's constant at given temperature and pressure (mol/kg/atm)
    real KH(real TK) {
    // using version from Plummer and Busenberg (1982)
    real kh = 10^(108.3865 + 0.01985076 * TK - 6919.53 / TK -
    40.45154 * log10(TK) + 669365 / TK^2);
    return kh;
    }
    
    // Apparent, or stoichiometric, solubility constant for calcite
    // Mucci 1983
    // Ca2++ + CO3-- = CaCO3 (s)
    real Ksp(real TK, real S) {
    // this is the -log(thermodynamic solubility constant) as a function of temp.
    // in distilled water according to Plummer and Busenberg 1982
    real pKsp_0 = -171.9065 - 0.077993 * TK + 2839.319 / TK + 71.595 * log10(TK);
    // These are how salinity affects the constant according to
    // Mucci 1983, Table 7
    real b = +(-0.77712 + 0.0028426 * TK + 178.34 / TK) * sqrt(S);
    real c = -0.07711 * S + 0.0041249 * S^1.5;
    real log10Kspc = pKsp_0 + b + c;
    real Kspc = 10^(log10Kspc);
    return Kspc;
    }
    
    // carbonate equilibrium function
    vector carb(real TK, real AT, real cond, real TC) {
    
    // Calculate all the parameters
    real gamma = act(TK, cond); // -log(monovalent act. coef)
    real I = 1.6 * 10^-5 * cond; // ionic strength, cond in uS/cm
    real S = 53.974 * I; //salinity from ionic strength, estimated
    real aH = 10 ^ -gamma; // activity coefficient for H+
    real aOH = 10 ^ -gamma; // activity coefficient for OH-
    real aHCO3 = 10 ^ -gamma; // activity coefficient for HCO3-
    real aCO3 = 10 ^ (4 * -gamma); // activity coefficient for CO32-
    real KWa = Kw(TK) / (aH * aOH); // apparent dissociation coefficient
    real K1a = K1(TK) / (aH * aHCO3); // apparent dissociation coefficient
    real K2a = K2(TK) / (aH * aCO3 / aHCO3); // apparent dissociation coefficient
    
    // Calculate pH
    // Simple iterative scheme from Park 1969 in L&O
    // Iterate for H and CA by repeated solution
    real pH = 8; // initial guess     
    real H = 10 ^ (-pH);  
    real delH = H;   
    real tol = 1.e-15; // tolerance for difference in iterations
    int iter = 0;
    
    while (delH > tol) {  // iterate until H converges
    real H_old = H;  // previous H
    // solve for carbonate alkalinity from TA
    real CA = AT / 1000;
    
    // solve quadratic for H
    real a = CA;
    real b = K1a * (CA - TC / 1000);
    real c = K1a * K2a * (CA - 2 * TC / 1000);
    H = (-b + sqrt(b^2 - 4 * a * c) ) / (2 * a);
    
    // How different is new estimate from previous one?
    delH = abs(H - H_old);
    iter = iter + 1;
    
    }
    pH = -log10(H);
    real CT = TC / 1000;
    
    // Solve the carbonate system, from Stumm and Morgan, (mol/kg)
    real alpha1 = (H * K1a) / (H^2 + K1a * H + K1a * K2a);  // HCO3 ionization fraction
    real alpha2 = (K1a * K2a) / (H^2 + K1a * H + K1a * K2a); // CO3 ionization fraction
    real CO2 = CT * (H ^ 2) / (H ^ 2 + K1a * H + K1a * K2a);
    real HCO3 = CO2 * K1a / H;
    real CO3 = HCO3 * K2a / H;
    
    // Get concentrations in mol/m3
    CO2 = CO2 * 1000;
    HCO3 = HCO3 * 1000;
    CO3 = CO3 * 1000;
    
    vector[4] sys;
    sys[1] = pH;
    sys[2] = CO2;
    sys[3] = HCO3;
    sys[4] = CO3;
    // Return the carbonate system
    return sys;
    }
      // Function that calculates the derivatives of the ODE
    vector carbODE(real t, vector y, vector theta, real par, real temp, vector dvals) {
        
    // Unpack states
    real O2  = y[1];
    real DIC = y[2];
    real ALK = y[3];
    
    // Unpack parameters
    real qL_w     = theta[1];
    real ALK_l    = theta[2];
    real DIC_l    = theta[3];
    real O2_l     = theta[4];
    real gpp_mean = theta[5];
    real er_mean  = theta[6];
    real k600     = theta[7];
    real PQ       = theta[8];
    real RQ       = theta[9];
    real z        = dvals[1];
    real cond     = dvals[2];
    real meanpar  = dvals[3];
   
    // Calculate Schmidt numbers for CO2 and O2 (Rik Wanninkhof L&O methods 2014)
    real Sc_CO2 = 1923.6 - 125.06 * temp + 4.3773 * temp^2 - 0.085681 * temp^3 + 0.00070284 * temp^4;
    real Sc_O2  = 1745.1 - 124.34 * temp + 4.8055 * temp^2 - 0.10115 * temp^3 + 0.00086842 * temp^4;
    
    // Calculate gas exchange coefficients for CO2 and O2 (m/hr)
    real k_CO2 = k600 * (600 / Sc_CO2)^0.5;
    real k_O2  = k600 * (600 / Sc_O2)^0.5;
    
    // // Calculate gas exchange velocity (m/h)
    // real k_CO2 = K_CO2 * z / (24);
    // real k_O2  = K_O2 * z / (24);
    
    // Temperature in Kelvin
    real TK = temp + 273.15;
    
    // -log(monovalent act. coef)
    real gamma = act(TK, cond);
    
    // Calculate Henry's constant K0 for CO2 (mol/kg/atm); Weiss 1974 Mar. Chem.
    // using freshwater version from Plummer and Busenberg (1982)
    real K_henry = KH(TK);
    
    // Calculate saturation value for O2
    real O2_sat = o2_sat(temp) / 32; // (mol/m3)
    
    // Calculate saturation value for CO2
    real CO2_sat = 400 * 1e-6 * 1 * K_henry * 1000; // (mol/m3)
    
    // Calculate carbonate equilibrium
    // compute pH, HCO3, and CO2 from alkalinity and DIC
    vector[4] carb_eq = carb(TK, ALK, cond, DIC);
    
    // Get the equilibrium values of the carbonate system
    // pH    = carb_eq[1];
    real CO2   = carb_eq[2];
    // HCO3  = carb_eq[3];
    // CO3   = carb_eq[4];
    // Ca    = HCO3 / 2; //assuming 1:2 relationship in calcium carbonate system

    // Calculate fluxes
    // First calculate GPP as linear function of mean PAR (mol O2/m2/hr)
    real gpp = (gpp_mean / 24 / 32) * (par / meanpar);
    
    // Ecosystem respiration (mol O2/m2/hr)
    real er = (er_mean / 24) / 32;
    
    // Gas exchange (mol/m2/del_t) from river perspective (+ is into river)
    real fO2  = k_O2 * (O2_sat - O2);
    real fCO2 = k_CO2 * (CO2_sat - CO2);
    
    // exchange with lateral inflow  (mol/m2/hr)
    real LO2  = qL_w * (O2_l - O2);
    real LDIC = qL_w * (DIC_l - DIC);
    real LALK = qL_w * (ALK_l - ALK);
    
    // Calculate derivatives (just took out z term)
    vector[3] dydt; // Number of state variables
    dydt[1] = (gpp - er + fO2 + LO2) / z; // Derivative of O2
    dydt[2] = (-gpp * (1 / PQ) + er * RQ + fCO2 + LDIC) / z; // Derivative of DIC
    dydt[3] = LALK / z; // Derivative of ALK
    
    return dydt;
    }

    // Function that performs the trapezoid rule
    vector trapezoid_rule(real t0, real t1, vector y0, vector theta, real par, real temp, vector dvals, int n_steps) {
        real dt = (t1 - t0) / n_steps;
        vector[3] y = y0;
        for (i in 1:n_steps) {
            vector[3] k1 = carbODE(t0 + (i - 1) * dt, y, theta, par, temp, dvals);
            vector[3] k2 = carbODE(t0 + i * dt, y + dt * k1, theta, par, temp, dvals);
            y += 0.5 * dt * (k1 + k2);
        }
        return y;
    }
}

data {
    int<lower=1> n_days;  // Number of days
    int<lower=1> n_obs;  // Number of observations per day
    int<lower=1> n_steps;  // Number of steps for trapezoid rule
    real<lower=0> dt;  // Time step for observations
    vector[n_days] t0;  // Start time for each day
    vector[n_days] t1;  // End time for each day
    array[n_obs, n_days] real O2_obs;  // Observed O2 data
    array[n_obs, n_days] real CO2_obs;  // Observed CO2 data
    array[n_obs, n_days] real par;  // Observed PAR
    array[n_obs, n_days] real temp;  // Observed temperature
    array[n_days] vector[3] dvals; // Daily water depth, conductivity, and mean PAR
    array[n_days] vector[3] y0_init;  // Initial states for each day
}

parameters {
    // Parameters to estimate with their bounds
    vector<lower=0>[n_days] qL_w;
    vector<lower=0>[n_days] ALK_l;
    vector<lower=0>[n_days] DIC_l;
    vector<lower=0>[n_days] O2_l;
    vector<lower=0>[n_days] gpp_mean;
    vector<lower=gpp_mean/2>[n_days] er_mean; //ER is at least half of GPP
    vector<lower=0>[n_days] k600;
    vector<lower=0.5, upper=2>[n_days] PQ;
    vector<lower=0.5, upper=2>[n_days] RQ;
    real<lower=0> sigma_y;  // Observation error standard deviation
    vector<lower=0>[9] mu_theta;
}

transformed parameters {
    matrix[n_obs, n_days] O2_hat; // we always have O2 observations to test against
    matrix[n_obs, n_days] CO2_hat; // and always CO2 observations to test against
    matrix[n_obs, n_days] pH_hat; // sometimes we have pH observations to test against
    vector[9] theta;
    for (d in 1:n_days) {
        theta[1] = qL_w[d];
        theta[2] = ALK_l[d];
        theta[3] = DIC_l[d];
        theta[4] = O2_l[d];
        theta[5] = gpp_mean[d];
        theta[6] = er_mean[d];
        theta[7] = k600[d];
        theta[8] = PQ[d];
        theta[9] = RQ[d];
        for (i in 1:n_obs) {
            real t = t0[d] + (i - 1) * dt;
            vector[3] y_hat = trapezoid_rule(t, t + dt, y0_init[d], theta, par[i, d], temp[i, d], dvals[d], n_steps);
            O2_hat[i, d] = y_hat[1];
            // Compute pH and CO2 using the carb function for each day and time
            vector[4] carb_hat = carb(temp[i, d] + 273.15, y_hat[3], dvals[d, 2], y_hat[2]);
            pH_hat[i, d] = carb_hat[1];
            CO2_hat[i, d] = carb_hat[2];
        }
    }
}

model {
   // Hyperprior distributions of parameters across days
    mu_theta[1] ~ normal(0.05, 0.02); // qL_w (m2/hr), covers a good range
    mu_theta[2] ~ normal(0.5, 0.3); // ALK_l (mol/m3), covers a good range
    mu_theta[3] ~ normal(1, 0.3); // DIC_l (mol/m3), covers a good range
    mu_theta[4] ~ normal(0.1, 0.05); // O2_l (mol/m3), covers a good range
    mu_theta[5] ~ normal(8, 4); // gpp_mean (g/m2/d), covers a good range
    mu_theta[6] ~ normal(8, 4); // er_mean (g/m2/d), covers a good range
    mu_theta[7] ~ lognormal(3, 0.4); // K600 (m/hr), covers a good range
    mu_theta[8] ~ normal(1, 0.3); // PQ (-), covers a good range
    mu_theta[9] ~ normal(1, 0.3); // RQ (-), covers a good range

    qL_w ~ normal(mu_theta[1], 1);
    ALK_l ~ normal(mu_theta[2], 1);
    DIC_l ~ normal(mu_theta[3], 1);
    O2_l ~ normal(mu_theta[4], 1);
    gpp_mean ~ normal(mu_theta[5], 1);
    er_mean ~ normal(mu_theta[6], 1);
    k600 ~ normal(mu_theta[7], 1);
    PQ ~ normal(mu_theta[8], 1);
    RQ ~ normal(mu_theta[9], 1);

    // Priors for observation error
    sigma_y ~ cauchy(0, 1);

    // Likelihood for observed data
    for (d in 1:n_days) {
        for (i in 1:n_obs) {
            real t = t0[d] + (i - 1) * dt;
            vector[3] y_pred = trapezoid_rule(t, t + dt, y0_init[d], theta, par[i, d], temp[i, d], dvals[d], n_steps);
            O2_obs[i, d] ~ normal(y_pred[1], sigma_y);
            // Compute CO2 using the carb function for each day and time
            vector[4] carb_pred = carb(temp[i, d] + 273.15, y_pred[3], dvals[d, 2], y_pred[2]);
            CO2_obs[i, d] ~ normal(carb_pred[2], sigma_y);
        }
    }
}

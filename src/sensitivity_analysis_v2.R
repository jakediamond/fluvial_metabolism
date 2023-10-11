# -------------------------------------
# Author: 
# Purpose: parameter estimation for lateral inflows and concentrations
# Date: 
# -------------------------------------
# Load libraries
library(FME)

# Load model and associated functions
source(file.path("src", "carbonate_model.R"))
source(file.path("src", "initialize_model.R"))

# Make the paramter ranges to check sensitivity of
parRanges <- data.frame(min = c(0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0 , 0 , 1 , 0.0), 
                        max = c(0.5, 2.0, 0.3, 2.0, 2.0, 2.0, 20, 20, 10, 0.5))
rownames(parRanges) <- c("qL_w", "DIC_l", "O2_l", "ALK_l", "PQ", "RQ",
                         "gpp_mean", "er_mean", "K600", "CaCO3")
parRanges

# Need this form for the sensitivity analysis
solvemodel <- function(x, parset = names(x)) {
  pars[parset] <- x
  
  # Equilibrium based on pH and ALK from parameters
  # carb_ini <- carb(TK = pars["temp"] + 273.15, AT = pars["Alk"], 
  #                  pH = pars["pH"], cond = pars["cond"])
  carb_ini <- seacarb::carb(flag = 24,
                            var1 = pars["CO2_atm"], var2 = pars["Alk"] / rhow(pars["temp"]),
                            pHscale = "F", k1k2 = "m06",
                            S = pars["cond"] * 1.6E-5 * 54,
                            T = pars["temp"], warn = "n")
  
  # Initial O2 concentration and carbonate system (mol/m3)
  O2_ini   <- as.numeric(O2_sat(pars["temp"])) / 31.998  
  DIC_ini  <- as.numeric(carb_ini$DIC * rhow(pars["temp"]))
  ALK_ini  <- as.numeric(carb_ini$ALK * rhow(pars["temp"]))
  CALC_ini <- as.numeric(pars["CaCO3"]) 
  
  # Initial conditions vector
  yini <- c(O2 = O2_ini,
            DIC = DIC_ini,
            ALK = ALK_ini,
            CALC = CALC_ini)
  
  ## ode solves the model by integration ...
  return(as.data.frame(ode(y = yini, times = times, func = model,
                           parms = pars)))
}

# Then we estimate the sensitivity to one parameter, CaCO3 (parameter 10), varying 
# its values according to a regular grid (dist=grid). The effect of that on DIC and ALK
# is estimated. To do this, the model is run 50 times (num=50). The system.time 
# is printed (in seconds):
print(system.time(
  sR <- sensRange(func = solvemodel,
                  parms = pars,
                  dist = "grid",
                  sensvar = c("DIC", "ALK"),
                  parRange = parRanges[10, ], num = 50)
))
head(summary(sR))

summ.sR <- summary(sR)
par(mfrow=c(2, 2))
plot(summ.sR, xlab = "time, hour", ylab = "molC/m3",
     legpos = "topright", mfrow = NULL)
plot(summ.sR, xlab = "time, hour", ylab = "molC/m3", mfrow = NULL,
     quant = TRUE, col = c("lightblue", "darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 3, "Sensitivity to CaCO3", cex = 1.25)
par(mfrow = c(1, 1))

# Sensitivity ranges can also be estimated for a combination of parameters. 
# Here we use all parameters, and use a latin hypercube sampling algorithm.
parms_vary <- pars[c("PQ", "RQ", "ALK_l", "O2_l", "DIC_l", "QL_f",
                           "gpp_mean", "er_mean", "K600", "CaCO3")]
Sens2 <- summary(sensRange(func = solvemodel,
                           parms = parms_vary,
                           dist = "latin", sensvar = "DIC",
                           parRange = parRanges, num = 50))
plot(Sens2, main = "Sensitivity", xlab = "time, hour",ylab = "molC/m3")



# Local sensitivity he effect of a parameter value in a very small region near 
# its nominal value is estimated. The methods are based on Brun et al. (2001)
Sns <- sensFun(func = solvemodel, parms = parms_vary,
               sensvar = "DIC", varscale = 1)
head(Sns)
plot(Sns)
# Based on the sensitivity functions, several summaries are generated, 
# which allow to rank the parameters based on their influence.
summary(Sns)
# Looks like HCO3_l, CO2_l, and RQ are the most sensitive
cor(Sns[ ,-(1:2)])
pairs(Sns)


# # Monte carlo runs
# SF <- function (pars) {
#   out <- solvmodel(pars)
#   return(out[nrow(out), 2:3])
#   }
# CRL <- modCRL(func = SF, parms = parms_data, parRange = parRanges[1,])
# plot(CRL)

# Based on the sensitivity functions of model variables to selection of parameters, function
# collin calculates the collinearity or identifiability of sets of parameters.
Coll <- collin(Sns)
Coll
Coll [Coll[,"collinearity"] < 20 & Coll[ ,"N"] == 4, ]


# # Monte Carlo methods can also be used to see how parameter uncertainties propagate, i.e. to
# # derive the distribution of output variables as a function of parameter distribution.
# # Here the effect of the parameters gmax and eff on final bacterial concentration is assessed. The
# # parameter values are generated according to a multi-normal distribution; they are positively
# # correlated (with a correlation = 0.63).
# CRL2 <- modCRL(func = SF, parms = pars, parMean = c(gmax = 0.5, eff = 0.7),
#                parCovar = matrix(nr = 2, data = c(0.02, 0.02, 0.02, 0.05)),
#                dist = "norm", sensvar = "DIC", num = 150)


# model fit ---------------------------------------------------------------
# Data to compare
Dat <- pluck(mod_data, 7, 11) %>%
  select(time_hr, O2, ALK = Alk, DIC) %>%
  mutate(time = time_hr / del_t) %>%
  select(-time_hr) %>%
  select(time, O2, DIC, ALK) %>% 
  filter(time <= sim_time) %>%
  as.data.frame()

# FME function modCost estimates the “model cost”, which the sum of (weighted) squared
# residuals of the model versus the data. This function is central to parameter identifiability
# analysis, model fitting or running a Markov chain Monte Carlo.

# We first define an objective function that returns the residuals of the model versus the data,
# as estimated by modcost. Input to the function are the current values of the parameters that
# need to be finetuned and their names (or position in par).
Objective <- function(x, parset = names(x)) {
  # parms_data[parset] <- x
  # tout <- seq(0, 50, by = 0.5)
  ## output times
  out <- solvemodel(x)
  ## Model cost
  return(modCost(obs = Dat, model = out, weight = "std"))
}


# obj_func <- function(params) {
#   predicted <- solvemodel(params)
#   return(sum((Dat[, 2:4] - predicted[, 2:4])^2))
# }
# First it is instructive to establish which parameters can be identified.
# We assess that by means of the function collin, selecting only the output
# variables at the instances when there is an observation.
Coll <- collin(sF <- sensFun(func = Objective, parms = parms_vary, 
                             varscale = 1))
Coll

# The larger the collinearity value, the less identifiable the parameter based on the data. In
# general a collinearity value less than about 20 is ”identifiable”. Below we plot the collinarity
# as a function of the number of parameters selected. We add a line at the height of 20, the
# critical value:
plot(Coll, log = "y")
abline(h = 20, col = "red")
Coll [Coll[,"collinearity"] < 20 & Coll[ ,"N"] == 6, ]
#don't really care about gpp_mean or K600 or er_mean or RQ

# We now use function modFit to locate the minimum. It includes several fitting procedures;
# the default one is the Levenberg-Marquardt algorithm.
# In the following example, parameters are constrained to be > 0
## start values for the parameters
plower    <- c(QL_f = 0, O2_l = 0.01, CO2_l = 0.1, HCO3_l = 1.5, 
               PQ = 0.8, CaCO3 = 0, er_mean = 8)
pupper    <- c(QL_f = 0.5,O2_l = 0.3,  CO2_l = 1, HCO3_l = 4, 
               PQ = 2, CaCO3 = 1, er_mean = 18)
print(system.time(Fit <- modFit(p = c(QL_f = 0.2, CO2_l = 0.2, O2_l = 0.3,
                                      HCO3_l = 2, PQ = 1, CaCO3 = 0.5, er_mean = 14),
                                f = Objective, upper=pupper, lower=plower,
                                # method = "Nelder-Mead")))
                                control = list(nprint = 1))))

pars
# Parameters to test
parms_test <- c(CO2_l = 0.015, HCO3_l = 2, PQ = 1.6, RQ = 2.7, K600 = 2.4, CaCO3 = 0.1)
# Lower bounds
plower    <- c(CO2_l = 0.01, HCO3_l = 1.5, PQ = 0.8, RQ = 0.8, K600 = 1, CaCO3 = 0.01)
# Upper bounds
pupper    <- c(CO2_l = 0.5, HCO3_l = 2.5, PQ = 2, RQ = 3, K600 = 4, CaCO3 = 0.5)

print(system.time(Fit4 <- modFit(p = parms_test,
                                f = Objective, upper=pupper, lower=plower,
                                control = list(nprint = 1))))
summary(Fit4)
parms_fit4 <- pars
parms_fit4[c("CO2_l", "HCO3_l", "PQ", "RQ", "K600", "CaCO3")] <- Fit4$par
out_fit4 <- solvemodel(parms_fit4)

plot(Dat$time, Dat$DIC)
lines(out_fit4$DIC, col = "red")
plot(Dat$time, Dat$ALK)
lines(out_fit4$ALK, col = "red")
plot(Dat$time, Dat$O2)
lines(out_fit4$O2, col = "red")

plot(Dat$DIC, Dat$O2)
points(out_fit4$ALK, out_fit4$O2, col = "red")
plot(out_fit3)
lm(O2~DIC, data = out_fit3)
lm(O2~DIC, data = Dat)
lm(O2~ALK, data = out_fit3)
lm(O2~ALK, data = Dat)

fit3_bfgs <- modFit(f = Objective, p = parms_test,
                              method="L-BFGS-B",
              upper = pupper, lower = plower, control = list(trace = 2,
                                                             REPORT = 1,
                                                             parscale = parms_test))
summary(fit3_bfgs)
parms_fit3_bfgs <- pars
parms_fit3_bfgs[c("CO2_l", "HCO3_l", "PQ", "RQ", "K600")] <- fit3_bfgs$par
out_fit3_bfgs <- solvemodel(parms_fit3_bfgs)

plot(Dat$time, Dat$DIC)
lines(out_fit3_bfgs$DIC, col = "red")
plot(Dat$time, Dat$ALK)
lines(out_fit3_bfgs$ALK, col = "red")
plot(Dat$time, Dat$O2)
lines(out_fit3_bfgs$O2, col = "red")

plot(out_fit3_bfgs)
lm(O2~DIC, data = out_fit3_bfgs)
lm(O2~DIC, data = Dat)
lm(O2~ALK, data = out_fit3_bfgs)
lm(O2~ALK, data = Dat)







parms_fit_og <- pars
parms_fit_og[c("QL_f", "CO2_l", "O2_l", "HCO3_l", "PQ", "CaCO3")] <- 
  c(QL_f = 0.25, CO2_l = 0.014, O2_l = 0.193, HCO3_l = 2.04, PQ = 1.28, CaCO3 = 0.025)

# init <- solvemodel(pars)
parms_fit <- pars
parms_fit[c("QL_f", "CO2_l", "O2_l", "HCO3_l", "PQ", "CaCO3", "er_mean")] <- Fit$par
parms_fit
out_fit <- solvemodel(parms_fit)
out_fit_og <- solvemodel(parms_fit_og)
# Cost <- modCost(obs = Dat, model = out_fit)
# Cost
# plot(Cost$residuals)
# plot(out_fit)

# plot(out_fit)
# plot(Cost, xlab = "time", ylab = "", main = "residuals")

plot(Dat$time, Dat$DIC)
lines(out_fit_og$DIC, col = "red")
plot(Dat$time, Dat$ALK)
plot(out_fit_og$ALK, col = "red")
plot(Dat$time, Dat$O2)
lines(out_fit$O2, col = "red")
# plot(init, out_fit, col = c("green", "blue"), lwd = 2, lty = 1, 
#      obs = Dat, obspar = list(col = "red", pch = 16, cex = 2), 
#      main = "model")

legend("right", c("initial", "fitted"), col = c("green", "blue"), lwd = 2)


parms_fit2 <- pars
parms_fit2[c("QL_f", "CO2_l", "O2_l", "HCO3_l", "PQ", "CaCO3", "RQ", "K600")] <- 
  c(QL_f = 0.25, CO2_l = 0.001, O2_l = 0.25, HCO3_l = 2, PQ = 1.28, 
    CaCO3 = 0.00, RQ = 2.5, K600 = 2.9)

out_fit2 <- solvemodel(parms_fit2)
# Cost <- modCost(obs = Dat, model = out_fit)
# Cost
# plot(Cost$residuals)
# plot(out_fit)

# plot(out_fit)
# plot(Cost, xlab = "time", ylab = "", main = "residuals")

plot(Dat$time, Dat$DIC)
lines(out_fit2$DIC, col = "red")
plot(Dat$time, Dat$ALK)
lines(out_fit2$ALK, col = "red")
plot(out_fit2$O2, col = "red")
points(Dat$time, Dat$O2)
plot(out_fit2$ALK, col = "red")
points(Dat$time, Dat$ALK)





## fit the model; nprint = 1 shows intermediate results
fit <- modFit(f = modelCost, p = c(QL_f =0.1, CO2_l = 0.1, O2_l = 0.01, 
                                   PQ = 1, RQ = 1),
              upper=pupper, lower=plower, control = list(nprint = 1))
summary(fit)

## graphical result
out2 <- ode(y = y, parms = startpars, times = times, func = derivs)
out3 <- ode(y = y, parms = fit$par, times = times, func = derivs)
plot(out, out2, out3, obs = yobs)

legend("topleft", legend=c("original", "startpars", "fitted"),
       col = 1:3, lty = 1:3)






# Parameters and inits for the model based on Loire
# This is a list of datasets to compare from
mod_data <- readRDS(file.path("results", "model_comparison_data.RDS"))

# This is the parameter set to check for the Loire at Dampierre 
parms_data <- unlist(pluck(mod_data, 12, 11))[c("PQ", "RQ", "CO2_l", 
                                                "O2_l", "HCO3_l", "QL_f")]
# The light signal
light_signal <- light_f(pluck(mod_data, 2, 11), pluck(mod_data,5, 11), pluck(mod_data, 6,11))

temp_signal <- temp_f(pluck(mod_data, 2, 11), pluck(mod_data,5, 11), pluck(mod_data, 6,11))
ini_data <- unlist(pluck(mod_data, 13, 11))
ini_data <- c(ini_data, CALC = 0.5)

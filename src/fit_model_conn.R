
library(tidyverse)
library(deSolve)
library(FME)


# Load the model and global parameters and initial conditions
source(file.path("src", "carbonate_model.R"))
source(file.path("src", "initialize_model.R"))
source(file.path("src", "compare_functions.R"))

df_data <- read_fun(file.path("data", "connecticut", "connecticut_clean.xlsx")) |>
  filter(site == "thom") |>
  mutate(datestart = "2020-09-07", dateend = "2020-09-07") |> #3/10 before
  mutate(
    across(datestart, as.Date),
    across(dateend, as.Date))

# Get all the model data necessary
df_moddat <- df_data |>
  mutate(# time step
    delt = map(timeseries, delt_f),
    # data to compare
    dat_comp = pmap(list(timeseries, datestart, dateend, delt), ts_f),
    # initial conditions of data
    inits = map2(timeseries, datestart, init_f),
    # daily metabolism data
    met = NA, #map2(daily, datestart, met_f),
    # light forcing function based on measurements
    lightfun = pmap(list(timeseries, datestart, dateend, delt), light_f),
    # temperature forcing function based on measurements
    tempfun = pmap(list(timeseries, datestart, dateend, delt), temp_f),
    # update model parameters with data
    mod_parms = pmap(list(inits, parms, met), update_parms),
    # update model initial conditions
    mod_inits = map2(inits, mod_parms, update_inits))

parms_data <- unlist(pluck(df_moddat, 13, 1))#[c("RQ", "PQ", "DIC_l",
#  "ALK_l", "qL")] #"PQ", "O2_l",

parms_data
colnames(df_moddat)
light_signal <- pluck(df_moddat, 11, 1)

temp_signal <- pluck(df_moddat, 12, 1)
temp_signal(1)
ini_data <- unlist(pluck(df_moddat, 14, 1))[c("O2", "DIC", "ALK")]

?seacarb::carb
seacarb::carb(flag = 5, var1 = 5.982E-2/1000, var2 = 0.807/1000, S=0.2, T = 23.51)

ini_data<- c(O2 = 0.235, DIC = 0.807, ALK = 0.749)
# ini_data <- c(O2 = 0.32, DIC = 1.82, ALK = 1.76)
datcomp <- pluck(df_moddat, 8, 1) |>
  select(time = time_hr,
         # O2,
         CO2,
         # pH,
         DIC
  ) |>
  as.data.frame()


plot(datcomp)

del_t    <- unlist(pluck(df_moddat, "delt"))                # time step (h)
days     <- round((max(datcomp$time) - min(datcomp$time)) / 24)                 # number of days to simulate
sim_time <- (24 * days ) - 1 *del_t # simulation time (h)
times    <- seq(0, sim_time, # sequence of in-model times to simulate
                by = del_t)  
# ini_data2 <- c(ini_data, "CALC" = parms_data[["CaCO3"]])
# ini_data2

# parms_data
# We first define an objective function that returns the residuals of the model versus the data,
# as estimated by modcost. Input to the function are the current values of the parameters that
# need to be fine tuned and their names (or position in par).
cost <- function(p, p_dat = parms_data, i_dat = ini_data, t = times, dat = datcomp,
                 mod = model, ls = light_signal, ts = temp_signal) {
  
  # The names of the parameters to be changed
  p_dat[names(p)] <- p
  i_dat <- c(i_dat, "CALC" = p_dat[["CaCO3"]])
  
  # The output of the model
  out <- ode(i_dat, t, mod, p_dat, 
             light_forcing = FALSE, temp_forcing = TRUE)
  
  # return(modCost(obs = dat, model = out))
  # sum((scale(dat[-1]) - scale(out1[-1]))^2)
  return(modCost(out, dat))
}

parms_test <- c(K600 = 4, DIC_l = 2.3, ALK_l = 2.3, O2_l = 0.2, qL = 4.6E-5, z = 3,
                gpp_mean = 15, er_mean = 15)
plower    <- c(K600 = 0.5, DIC_l = 0.5, ALK_l = 0.5, O2_l = 0, qL = 1E-5, z = 2.9,
               gpp_mean = 6, er_mean = 6)
pupper    <- c(K600 = 8, DIC_l = 3, ALK_l = 3, O2_l = 0.4, qL = 1E-3, z = 4,
               gpp_mean = 35, er_mean = 30)
# Limited-memory Broyden–Fletcher–Goldfarb–Shanno algorithm with box constraints
# (L-BFGS-B) is a quasi-Newton method optimization algorithm. It uses a limited
# amount of computer memory in a linear fashion and works well for problems
# with many variables.
fit2 <- modFit(f = cost, p = parms_test, method="L-BFGS-B",
               upper = pupper, lower = plower, control = list(trace = 1,
                                                              REPORT = 1,
                                                              parscale = parms_test))
summary(fit2)
parms_fit2 <- replace(parms_data, names(coef(fit2)), coef(fit2))
i_dat_fit2 <- c(ini_data, "CALC" = parms_data[["CaCO3"]])
fitout2 <- ode(i_dat_fit2, times, model, parms_fit2, 
               light_forcing = FALSE, temp_forcing = TRUE)
parms_fit2
# out2 <- ode(ini_data, times, model_to_solve2, coef(fit))
plot(fitout2, obs=datcomp, obspar=list(pch=16, col="red"))
 # x <- ode_to_df(fitout2)
# par <- par_fun(0:23, day = 251, latitude = 41.99)
# par
# plot_fun(ode_to_df(fitout2))
var0 <- fit2$var_ms_unweighted
cov0 <- summary(fit2)$cov.scaled * 2.4^2/5

MCMC <- modMCMC(f = cost, p = fit2$par,upper = pupper, lower = plower,
                niter = 5000)#, jump = 8,
                # var0 = var0, wvar0 = 0.1, updatecov = 50)

summary(MCMC)
plot(MCMC)

# These are the results from a near fit
pfit2    <- c(K600 = 0.95, DIC_l = 1.97, ALK_l = 0.88, O2_l = 0.31,
              gpp_mean = 10.2, er_mean = 5, PQ = 1.5, RQ = 1.5)

# These are the results from a near fit
pfit    <- c(K600 = 1.92, DIC_l = 2.11, ALK_l = 0.66, 
             gpp_mean = 14.4, er_mean = 5, PQ = 2, RQ = 1.89)
# Compare pco2 measured to estimated just from equilibria from clark flork
df_cf <- readxl::read_xlsx(file.path("data", "montana", "montana_clean.xlsx"))
# str(df_cf)
df_cf <- df_cf %>%
  mutate(pCO2_est = carb(temp + 273.15,
                         ALK,
                         pH,
                         cond)$pCO2)
# pCO2_est2 = carb2(temp + 273.15,
#                  ALK,
#                  pH,
#                  cond)$pCO2,
# pCO2_cor = carb2(temp + 273.15,
#                  ALK,
#                  pH,
#                  cond)$pCO2_cor)

ggplot(data = df_cf,
       aes(x = pCO2,
           y = pCO2_est)) +
  geom_point() +
  ggpubr::stat_regline_equation() +
  theme_classic() +
  labs(x = expression(pCO[2]~"measu. (ppm)"),
       y = expression(pCO[2]~"equil. (ppm)"))

ggplot(data = df_cf,
       aes(x = datetime)) +
  geom_point(aes(y = pCO2)) +
  geom_point(aes(y = pCO2_est), color = "red") +
  geom_line(aes(y = pCO2_est), color = "red") +
  scale_x_datetime(date_breaks = "1 day") +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "black")) +
  # ggpubr::stat_regline_equation() +
  
  labs(x = "",
       y = expression(pCO[2]~"(ppm)"))


# 
# ggplot(data = df_cf,
#        aes(x = pCO2,
#            y = pCO2_cor)) +
#   geom_point() +
#   ggpubr::stat_regline_equation() +
#   theme_classic() +
#   labs(x = expression(pCO[2]~"measu. (ppm)"),
#        y = expression(pCO[2]~"equil. (ppm)"))

parms_test <- c(RQ = 1.3, PQ = 1, Alk = 1.5)
## Now plot original and fitted models and data
out1 <- ode(ini_data, times, model_to_solve, parms_test)
out2 <- ode(ini_data, times, model_to_solve2, coef(fit))
plot(out2, obs=dat, obspar=list(pch=16, col="red"))


plower    <- c(RQ = 0.8, PQ = 0.8, Alk = 1)
pupper    <- c(RQ = 1.3, PQ = 1.3, Alk = 3)
cl <- makeCluster(4)     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster
clusterEvalQ(cl, library("deSolve"))
clusterCall(cl, light_signal)

fit_parallel <- optimParallel(f = cost, p = parms_test, method="L-BFGS-B",
                              upper = pupper, lower = plower, 
                              p_dat = parms_data,
                              i_dat = ini_data,
                              t = times,
                              mod = model_to_solve,
                              ls = light_signal,
                              ts = temp_signal,
                              control = list(trace = 2,
                                             REPORT = 1,
                                             parscale = parms_test),
                              parallel=list(loginfo=TRUE))
summary(fit_parallel)






# 
# parms_test <- c(K600 = 1, O2_l = 0.1, DIC_l = 1, ALK_l = 1, qL = 0.000329,
#                 gpp_mean = 6, er_mean = 12, PQ = 1, RQ = 1)
# 
# print(system.time(
#   fit <- modFit(f = cost, p = parms_test)
# ))
# summary(fit)
# # plot(fit)
# # ini_data_fit
# parms_fit <- replace(parms_data, names(coef(fit)), coef(fit))
# i_dat_fit <- c(ini_data, "CALC" = parms_data[["CaCO3"]])
# # ini_data_fit <- c(ini_data, "CALC" = coef(fit)[["CaCO3"]])
# fitout <- ode(i_dat_fit, times, model, parms_fit, 
#               light_forcing = TRUE, temp_forcing = TRUE)
# 
# plot(fitout, obs=datcomp, obspar=list(pch=16, col="red"))
# 
# output <- ode_to_df(fitout)
# 
# plot_fun(output, "CALC")
# 
# ggplot() +
#   geom_line(data = output,
#             aes(x = time,
#                 y = SI))
# # geom_point(data = rename(x, time = time_hr),
# #            aes(x = time,
# #                y = O2))
# 
# parms_test_o2 <- c(K600 = 4, O2_l = 0.2, qL = 4.6E-5,
#                    gpp_mean = 15, er_mean = 15)
# plower_o2 <- c(K600 = 1, O2_l = 0.01, qL = 1E-5,
#                gpp_mean = 5, er_mean = 5)
# pupper_o2 <- c(K600 = 8, O2_l = 0.2, qL = 1E-4,
#                gpp_mean = 30, er_mean = 20)
# fito2 <- modFit(f = cost, p = parms_test_o2, method="L-BFGS-B",
#                 upper = pupper_o2, lower = plower_o2, control = list(trace = 1,
#                                                                      REPORT = 1,
#                                                                      parscale = parms_test_o2))
# summary(fito2)
# parms_fito2 <- replace(parms_data, names(coef(fito2)), coef(fito2))
# i_dat_fito2 <- c(ini_data, "CALC" = parms_data[["CaCO3"]])
# fitouto2 <- ode(i_dat_fito2, times, model, parms_fito2, 
#                 light_forcing = FALSE, temp_forcing = TRUE)
# parms_fito2
# # out2 <- ode(ini_data, times, model_to_solve2, coef(fit))
# plot(fitouto2, obs=datcomp, obspar=list(pch=16, col="red"))
# 


# costs <- modCost(out1, dat, weight = "std")
# costs
# sum(costs)
# sum((dat[-1] - out1[-1])^2)
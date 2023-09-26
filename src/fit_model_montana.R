
library(tidyverse)
library(deSolve)
library(FME)


# Load the model and global parameters and initial conditions
source(file.path("src", "carbonate_model.R"))
source(file.path("src", "initialize_model.R"))
source(file.path("src", "compare_functions.R"))

df_data <- read_fun(file.path("data", "montana", "montana_clean.xlsx")) |>
  mutate(datestart = "2019-09-23", dateend = "2019-09-24") |>
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

# light_signal <- pluck(df_moddat, 11, 1)

temp_signal <- pluck(df_moddat, 12, 1)
ini_data <- unlist(pluck(df_moddat, 14, 1))[c("O2", "DIC", "ALK")]

datcomp <- pluck(df_moddat, 8, 1) |>
  select(time = time_hr,
         DIC,
         ALK,
         CO2) |>
  as.data.frame()

# ini_data2 <- c(ini_data, "CALC" = parms_data[["CaCO3"]])
# ini_data2

# parms_data
# We first define an objective function that returns the residuals of the model versus the data,
# as estimated by modcost. Input to the function are the current values of the parameters that
# need to be fine tuned and their names (or position in par).
cost <- function(p, p_dat = parms_data, i_dat = ini_data, t = times, 
                 mod = model, ls = light_signal, ts = temp_signal) {
  
  # # Forcing functions, need to be global
  # light_signal <<- ls
  # temp_signal <<- ts

  # The names of the parameters to be changed
  p_dat[names(p)] <- p
  i_dat <- c(i_dat, "CALC" = p_dat[["CaCO3"]])
  
  # The output of the model
  out <- ode(i_dat, t, mod, p_dat, 
             light_forcing = FALSE, temp_forcing = TRUE)
  
  return(modCost(obs = datcomp, model = out))
  # sum((scale(dat[-1]) - scale(out1[-1]))^2)
  # sum(modCost(out, dat, weight = "std")$var$SSR)
}

parms_test <- c(CaCO3 = 0.2, DIC_l = 3.2, ALK_l = 3.2, 
                gpp_mean = 5, er_mean = 5)

print(system.time(
fit <- modFit(f = cost, p = parms_test)
))
summary(fit)
plot(fit)

parms_fit <- replace(parms_data, names(coef(fit)), coef(fit))
ini_data_fit <- c(ini_data, "CALC" = coef(fit)[["CaCO3"]])
fitout <- ode(ini_data_fit, times, model, parms_fit, 
    light_forcing = FALSE, temp_forcing = TRUE)

# out2 <- ode(ini_data, times, model_to_solve2, coef(fit))
plot(fitout, obs=datcomp, obspar=list(pch=16, col="red"))

output <- ode_to_df(fitout)

plot_fun(output, "CALC")

ggplot() +
  geom_line(data = output,
            aes(x = time,
                y = SI))
  # geom_point(data = rename(x, time = time_hr),
  #            aes(x = time,
  #                y = O2))




# costs <- modCost(out1, dat, weight = "std")
# costs
# sum(costs)
# sum((dat[-1] - out1[-1])^2)
parms_test <- c(CaCO3 = 0.2, DIC_l = 3.2, ALK_l = 3.2, 
                gpp_mean = 5, er_mean = 5)#, qL = 3.5E-5)
plower    <- c(CaCO3 = 0.1, DIC_l = 2.8, ALK_l = 2.8, 
               gpp_mean = 2, er_mean = 2)#, qL = 2E-5)
pupper    <- c( CaCO3 = 0.5, DIC_l = 3.5, ALK_l = 3.5,
               gpp_mean = 8, er_mean = 8)#, qL = 5E-5)

fit2 <- modFit(f = cost, p = parms_test, method="L-BFGS-B",
              upper = pupper, lower = plower, control = list(trace = 1,
                                                             REPORT = 1,
                                                             parscale = parms_test))
summary(fit2)
parms_fit2 <- replace(parms_data, names(coef(fit2)), coef(fit2))
ini_data2 <- c(ini_data, "CALC" = coef(fit2)[["CaCO3"]])
fitout2 <- ode(ini_data2, times, model, parms_fit2, 
              light_forcing = FALSE, temp_forcing = TRUE)

# out2 <- ode(ini_data, times, model_to_solve2, coef(fit))
plot(fitout2, obs=datcomp, obspar=list(pch=16, col="red"))
plot(ode_to_df(fitout2)$SI)

plot_fun(ode_to_df(fitout2))


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
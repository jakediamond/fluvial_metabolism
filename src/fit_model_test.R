

# Load the model and global parameters and initial conditions
source(file.path("src", "carbonate_model.R"))
source(file.path("src", "initialize_model.R"))

df_data <- read_fun(file.path("data", "florida", "florida_clean.xlsx")) |>
  filter(site == "ICHE2700") |>
  mutate(datestart = "2019-07-09", dateend = "2019-07-10") |>
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
    met = map2(daily, datestart, met_f),
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
light_signal <- pluck(df_moddat, 11, 1)

temp_signal <- pluck(df_moddat, 12, 1)
ini_data <- unlist(pluck(df_moddat, 14, 1))[c("O2", "DIC", "ALK")]
ini_data <- c(ini_data, "CALC"=0)
datcomp <- pluck(df_moddat, 8, 1) |>
  select(time = time_hr,
         O2,
         CO2) |>
  as.data.frame()

parms_data
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
  
  # The output of the model
  out <- ode(i_dat, t, mod, p_dat, 
             light_forcing = TRUE, temp_forcing = TRUE)
  
  return(modCost(obs = datcomp, model = out))
  # sum((scale(dat[-1]) - scale(out1[-1]))^2)
  # sum(modCost(out, dat, weight = "std")$var$SSR)
}

parms_test <- c(RQ = 1.5, PQ = 1.2, DIC_l = 2, O2_l = 0.1, ALK_l = 1)

print(system.time(
fit <- modFit(f = cost, p = parms_test)
))
summary(fit)
plot(fit)
# costs <- modCost(out1, dat, weight = "std")
# costs
# sum(costs)
# sum((dat[-1] - out1[-1])^2)
p


fit <- modFit(f = cost, p = parms_test, method="L-BFGS-B",
              upper = pupper, lower = plower, control = list(trace = 1,
                                                             REPORT = 1,
                                                             parscale = parms_test))
summary(fit)
fit
plot(fit)
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
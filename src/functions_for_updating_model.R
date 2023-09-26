#
# Authors: Jake Diamond
# Purpose: Functions to allow rapid updating of model parameter sets 
# Date: 2023 April 1
# 
# function to get initial conditions based on date range
init_cond_fun <- function (timeseries, date_start) {
  filter(timeseries, date(datetime) == date_start) %>%
    head(1)
}

# function to get a timeseries from a date range
timeseries_fun <- function (timeseries, date_start, date_end) {
  filter(timeseries, date(datetime) >= date_start & date(datetime) <= date_end)
}

# function to replace x and y in the temp_signal() function based on date range
# the function should return a function in the global environment that can be 
# used in the model function. The function should calculate a variable called
# del_t that is the time step in hours. This time step is used to calculate
# the x values for the approxfun() function so they match the time step of the
# model
light_signal_fun <- function (timeseries, date_start, date_end) {
  # Get the time step in hours
  del_t <- timeseries$datetime[2] - timeseries$datetime[1]
  del_t <- as.numeric(del_t, units = "hours")
  # Get the x and y values for the approxfun() function
  x <- timeseries_fun(timeseries, date_start, date_end)$time_hr / del_t
  y <- timeseries_fun(timeseries, date_start, date_end)$light
  # Create the approxfun() function
  approxfun(x = x, 
            y = y, 
            method = "linear", rule = 2)
}

# Function to run the model with choices of parameters
update_parameters <- function(...){
  # Get the changes to parameters
  # arguments <- list(...)
  replace(parms,
          names(...),
          unlist(...))
}

# Function to run the model with choices of parameters
update_inis <- function(...){
  # Get the changes to parameters
  # arguments <- list(...)
  replace(yini,
          names(...),
          unlist(...))
}

mod_fun <- function (parameters, initial_conditions) {
  # Model with user specifications
  ode(y = initial_conditions,
      times = times,
      func = model,
      parms = parameters,
      temp_forcing = T,
      light_forcing = T)
}



# Function to run the model
mod_fun <- function (parameters, initial_conditions) {
  # Model with user specifications
  ode(y = unlist(initial_conditions),
      times = times,
      func = model,
      parms = unlist(parameters))
}

plot_fun_all <- function(data) {
  data <- tail(data, 24)
  p1 <- plot_fun(data, c("DIC", "exO2"), "biplot")
  p2 <- plot_fun(data, c("exCO2", "exO2"), "biplot")
  p <- (p1 + p2) + plot_layout(guides = "collect")
  # p3 <- plot_fun(day, c("ALK", "exO2"), "biplot")
  # p4 <- plot_fun(day, c("DIC", "ALK"), "biplot")
  # p <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")
}

lm_fun <- function(data) {
  day <- tail(data, 24) %>%
    mutate(dO2 = O2 - lag(O2),
           dCO2 = CO2 - lag(CO2),
           dDIC = DIC - lag(DIC),
           dALK = ALK - lag(ALK))
  # day <- filter(data, between(time_hr, 11, 17))
  o2_dic <- lm(O2 ~ DIC, data = day)
  o2_co2 <- lm(O2 ~ CO2, data = day)
  o2_alk <- lm(O2 ~ ALK, data = day)
  alk_dic <- lm(ALK ~ DIC, data = day)
  tibble(coef = c("int", "slope"),
         o2_dic = coef(o2_dic),
         o2_co2 = coef(o2_co2),
         o2_alk = coef(o2_alk),
         alk_dic = coef(alk_dic),
         calc = mean(day$calc),
         alk = mean(day$ALK),
         ql = mean(day$LALK)) %>%
    pivot_longer(cols = contains("_"))
}





# 
# # function to get the fractional hour of the day from a posixct object
# frac_hour <- function (datetime) {
#   hour(datetime) + minute(datetime) / 60 + second(datetime) / 3600
# }
# 
# 


# 
# # Stream temperature
# temp_signal <- approxfun(x = data$time_hr / del_t, 
#                          y = data$temp, 
#                          method = "linear", rule = 2)
# # Stream light
# light_signal <- approxfun(x = data$time_hr / del_t, 
#                           y = data$light, 
#                           method = "linear", rule = 2)

# 
# 
# # Update parameters function
# update_parms <- function(parm_trts, replace_vec = "pars") {
#   replace_vec <- get(replace_vec)
#   idx <- intersect(names(parm_trts), names(replace_vec))
#   # parm_trts[idx][is.na(parm_trts[idx])] <- y[which(is.na(dat[idx]))]
#   replace_vec[idx] <- parm_trts[idx]
#   replace_vec$ALK_l <- replace_vec$Alk
#   replace_vec$DIC_l <- replace_vec$Alk
#   replace_vec$Ca    <- replace_vec$Alk / 2
#   replace_vec
# }
# 
# # Update initial conditions functions
# update_inits <- function(params, replace_vec = "yini"){
#   replace_vec <- get(replace_vec)
#   # idx <- intersect(names(inis), names(replace_vec))
#   # replace_vec[idx] <- inis[idx]
#   temp <- as.numeric(params["temp"])
#   
#   # Equilibrium based on pH and ALK from parameters
#   carbeq <- seacarb::carb(flag = 24, 
#                           var1 = as.numeric(params["CO2_atm"]),
#                           var2 = as.numeric(params["Alk"]) / rhow(temp),
#                           T = temp,
#                           pHscale = "F",
#                           k1k2 = "m06",
#                           S = 1.5E-5*as.numeric(params["cond"])*54,
#                           warn = "n")
#   
#   # Stream DO concentration and carbonate system
#   replace_vec$O2  <- as.numeric(O2_sat(temp)) / 31.998
#   replace_vec$DIC <- as.numeric(carbeq$DIC * rhow(temp))# mol/m3
#   replace_vec$ALK <- as.numeric(params["Alk"]) # mol/m3
#   replace_vec$CALC <- as.numeric(params["CaCO3"]) # mol/m3
#   return(replace_vec)
# }
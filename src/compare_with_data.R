# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
# Load libraries
library(tidyverse)


# Choose a date range
datestart <- ymd("2019-05-07")
dateend   <- ymd("2019-05-08")
x = pluck(data, 2, 1)
init_cond_fun(x, ymd("2019-05-07"))

# function to get initial conditions based on date range
init_cond_fun <- function (timeseries, date_start) {
  filter(timeseries, date(datetime) == date_start) %>%
    head(1)
}

x <- data %>%
  drop_na() %>%
  mutate(inits = map(timeseries, ~init_cond_fun(., datestart)))


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

y = x %>%
  mutate(parameters = map(values, update_parameters),
         inits = map(inis, update_inis))

z <- y %>%
  mutate(mo = map2(parameters, inits, mod_fun))

# plot(z$mo[[2]])
# parms
# y$z
# y$values
# unlist(parameter_sets)
# names((x$values[[1]]))
# Create dataframe of treatments
parameter_sets <- data.frame(
  names   = c("loire_20220709", "ichetucknee_20190907"),
  latitude= c(47.8 , 30.0),
  doy     = c(190  , 190),
  z       = c(1    , 1),
  CO2_atm = c(418  , 400),
  p_atm   = c(1    , 1),
  temp    = c(20   , 22.1),
  pH      = c(8.39 , 7.24),
  Alk     = c(1.57 , 2.8), 
  Ca      = c(27.7 , 52),
  cond    = c(228  , 350),
  gpp_mean= c(1.3  , 1.33), 
  K600    = c(1.6  , 2.6),   
  er_mean = c(0.5  , 0.47)
)
27.7/1000/40.078

parms3 <-
  parms <- c(
    # geographic parameters
    latitude = 47.8,
    doy = 308, # day of the year
    
    # Hydraulic parameters
    z = 1, # (m) average depth
    
    # Physicochemical parameters
    CO2_atm = 391,   # (uatm) atmospheric pressure of CO2
    p_atm   = 1,     # (atm) barometric atmospheric pressure
    temp    = 21.5,    # (deg C) temperature of water
    pH      = 7.25,     # pH of the water
    Alk     = 2.9,   # (mol m-3)total measured alkalinity 
    Ca      = 1.3E-3,# (mol m-3) calcium concentration
    cond    = 333,   # (uS/cm) specific conductance
    
    # Reaction parameters
    gpp_mean = 9,  # (g O2/m2/d) daily mean GPP rate
    K600     = 2,   # (1/d) gas exchange coefficient at Schmidt number 600
    er_mean  = 9,  # (g O2/m2/d) daily ER rate
    PQ       = 1,   # photosynthetic quotient relating mols CO2 to produce mols O2
    RQ       = 1    # respiratory quotient relating mols O2 to produce mols CO2
  )

# Initial conditions vector
yini3 <- c(O2 = 0.14,
           DIC = 3.1,
           ALK = 2.93)

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
ode_out3 <- ode(y = yini3,
               times = times,
               func = model,
               parms = parms3,
               light_forcing = F,
               temp_forcing = F)

# Examine the model output ------------------------------------------------
# Get into good formats
output3 <- ode_to_df(ode_out3)
plot_fun(output3, "temp")












# Apply model to the treatments
mods <- trts %>%
  mutate(out = pmap(list(trts, func_mod)))

doini = 3.47

# Define the function that updates parameters based on user input
update_parameters <- function(parameters, ...) {
  # Update the specified parameters
  new_params <- list(...)
  for (param in names(new_params)) {
    if (param %in% names(parameters)) {
      parameters[[param]] <- new_params[[param]]
    } else {
      warning(paste0("Parameter '", param, "' not found in model."))
    }
  }
  
  # Return the updated parameters
  return(parameters)
}

# Define the initial parameter values
initial_params <- list(
  beta = 0.2,
  gamma = 0.1,
  mu = 0
)

# Define a series of parameter sets with names
parameter_sets <- list(
  set_names(list(pH = 8.39), "loire"),
  set_names(list(pH = 7.4), "ichetucknee")
)

# Run the model for each parameter set using purrr
results <- purrr::map(parameter_sets, function(params) {
  # Update the initial parameter values based on user input
  updated_params <- update_parameters(parms, !!!params)})
#   
#   # Define the initial state
#   initial_state <- c(S = 0.99, I = 0.01, R = 0)
#   
#   # Define the time vector
#   times <- seq(0, 100, by = 1)
#   
#   # Run the ODE model
#   output <- ode(y = initial_state, times = times, func = ode_model, parms = updated_params)
#   
#   # Return the output
#   return(output)
# })

























# Load data
df_loire <- readRDS(file.path("data", "loire", "all_hourly_data_complete.RDS"))

dat_comp <- filter(df_loire,
                   between(datetime, ymd_h(2022050900), ymd_h(2022051023))) %>%
  mutate(time_hr = row_number() - 1)

colnames(dat_comp)

ggplot() +
  geom_line(data = dat_comp,
            aes(x = time_hr / 24,
                y = CO2_uM)) +
  geom_line(data = output2,
            aes(x = time_hr / 24,
                y = CO2*1000),
            color = "red") +
  theme_bw()

plot_fun(output2, "GPP")



df_iche <- read_csv(file.path("data", "florida", "Ichetucknee", 
                              "ICHE2700_2019_15min_data2.csv")) %>%
  mutate(solartime = solar.time - hours(6))

dat_comp_ich <- filter(df_iche,
                       between(solartime, ymd_h(2019070900), ymd_h(2019071023))) %>%
  mutate(time_hr = (row_number() - 1) / 4)

ggplot() +
  geom_line(data = dat_comp_ich,
            aes(x = time_hr / 24,
                y = CO2_ppm/0.036/1000000)) +
  geom_line(data = output,
            aes(x = time_hr / 24,
                y = CO2),
            color = "red") +
  theme_bw()
colnames(dat_comp_ich)

35/12000

ceiling(25/24)
df_ich %>%
  filter(name %in% c("GPP", "ER")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(day = ceiling(time_hr / 24)) %>%
  group_by(day) %>%
  summarize(gpp = sum(GPP)*32,
            er = sum(ER) * 32)


plot(par(times * del_t))

henry(23+273.15)



output2 %>%
  mutate(day = ceiling(time_hr / 24)) %>%
  group_by(day) %>%
  summarize(gpp = sum(GPP)*32,
            er = sum(ER) * 32)

output2 %>%
  mutate(day = ceiling(time_hr / 24)) %>%
  group_by(day) %>%
  mutate(gpp2 = (10/32) * PAR/sumpar) %>%
  summarize(gpp = sum(GPP)*32,
            gpp3 = sum(gpp2) * 32,
            er = sum(ER) * 32)


















convert_dataframe <- function(df) {
  # Extract the column names of the first 9 columns
  col_names <- names(df)[2:14]
  
  # Combine the values of the first 9 columns into a list using pmap
  df_list <- pmap(df[, col_names], list) %>%
    map(unlist)
  
  # Create a new data frame with the name and list columns
  new_df <- tibble(
    name = df$name,
    values = df_list
  )
  
  return(new_df)
}
x = convert_dataframe(parameter_sets)













extract_columns(df_loire, ymd("20200702"))


extract_columns <- function(df, datetime, cols = c("O2", "pH")) {
  # Convert the datetime input to a column index
  datetime_col <- match(datetime, names(df))
  
  # Extract the columns based on the column names
  output_df <- df[, cols, drop = FALSE]
  
  # Check if any columns were not found
  not_found_cols <- setdiff(cols, names(df))
  
  # Generate warnings for any columns that were not found
  if (length(not_found_cols) > 0) {
    message("Warning: the following columns were not found in the data frame:")
    for (col in not_found_cols) {
      message(col)
    }
  }
  
  # Return the output data frame
  return(output_df)
}

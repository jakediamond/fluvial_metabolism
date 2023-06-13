# -------------------------------------
# Author: Jake Diamond
# Purpose: functions to prep measured data for use in the model
# Date: 2023-03-15
# -------------------------------------

# Function to load the data
read_fun <- function(file_path) {
  sheets <- set_names(readxl::excel_sheets(file_path))
  
  # Read each sheet into a data frame
  dfs <- map_dfr(sheets, ~ readxl::read_xlsx(file_path, sheet = .x, 
                                             guess_max = 2000) %>%
                   group_by(site) %>% nest(.key = .x) %>% ungroup())
  df <- dfs %>%
    pivot_longer(cols = -site) %>%
    drop_na() %>%
    pivot_wider()
}

# function to get delta t from data
delt_f <- function(data) {
  dt <- data$datetime[2] - data$datetime[1]
  dt <- as.numeric(dt, units = "hours")
  return(dt)
}

# function to get a timeseries from a date range
# also adds a column, "time_hr", which is continuous time in hours
# since the start of the timeseries
ts_f <- function (ts, date_start, date_end, dt) {
  dplyr::filter(ts, between(date(datetime), date_start, date_end)) %>%
    mutate(time_hr = (row_number() - 1) * dt)
}

# Generates a forcing light signal for the model based on a date range
light_f <- function (data, date_start, date_end, dt) {
  # Get the x and y values for the approxfun() function
  x <- ts_f(data, date_start, date_end, dt)$time_hr / dt
  y <- ts_f(data, date_start, date_end, dt)$light
  # Create the approxfun() function
  approxfun(x = x,
            y = y,
            method = "linear", rule = 2)
}

# Generates a forcing temp signal for the model based on a date range
temp_f <- function (data, date_start, date_end, dt) {
  # Get the x and y values for the approxfun() function
  x <- ts_f(data, date_start, date_end, dt)$time_hr / dt
  y <- ts_f(data, date_start, date_end, dt)$temp
  # Create the approxfun() function
  approxfun(x = x, 
            y = y, 
            method = "linear", rule = 2)
}

# Update parameters function
update_parms <- function(x, y, z, replace_vec = "pars") {
  replace_vec <- get(replace_vec)
  dat <- c(x, y, z)
  idx <- intersect(names(dat), names(replace_vec))
  dat[idx][is.na(dat[idx])] <- y[which(is.na(dat[idx]))]
  replace_vec[idx] <- dat[idx]
  replace_vec$doy <- yday(dat$date) 
  replace_vec
}

# Update initial conditions functions
update_inits <- function(inis, params, replace_vec = "yini"){
  replace_vec <- get(replace_vec)
  idx <- intersect(names(inis), names(replace_vec))
  replace_vec[idx] <- inis[idx]
  temp <- as.numeric(inis["temp"])
  # Equilibrium based on pH and ALK from parameters
  carbeq <- seacarb::carb(flag = 1,
                          var1 = as.numeric(inis["pH"]),
                          var2 = as.numeric(inis["CO2"]) / rhow(temp),
                          #flag = 8,
                          # var1 = as.numeric(params["pH"]),
                          # var2 = as.numeric(params["Alk"]) / 1000,
                          T = temp,
                          k1k2 = "m06",
                          S = 0,#params["cond"] * 1.6E-5 * 54,
                          warn = "n")
  
  # Stream DO concentration and carbonate system
  replace_vec$DIC <- carbeq$DIC * rhow(temp)   # mol/m3
  replace_vec$ALK <- carbeq$ALK * rhow(temp) #as.numeric(params["Alk"]) # mol/m3
  return(replace_vec)
}

# function to get initial conditions based on date range
init_f <- function (ts, date_start) {
  # First get the inital data from the timeseries if available
  inits <- filter(ts, date(datetime) == date_start) %>%
    head(1)
}

# function to get metabolism/daily conditions based on date range
met_f <- function (data, date_start) {
  dplyr::filter(data, date == date_start) %>%
    rename(gpp_mean = GPP, er_mean = ER) %>%
    mutate(er_mean = -er_mean)
}

# Function to run the model
mod_fun <- function (parameters, initial_conditions, driving_data, dt,
                     light, temp, 
                     l_force = FALSE, t_force = FALSE) {
  light_signal <<- light
  temp_signal <<- temp
  del_t <<- dt
  if("light" %in% colnames(driving_data)) {
    l_force = TRUE
  }
  if("temp" %in% colnames(driving_data)) {
    t_force = TRUE
  }
  times <- 1:nrow(driving_data)
  # Model with user specifications
  ode(y = unlist(initial_conditions),
      times = times,
      func = model,
      parms = unlist(parameters),
      temp_forcing = t_force,
      light_forcing = l_force)
}

# function to plot model and measurements
plot_comp_fun <- function (model_output, msmts, units = units_df) {

  bind_rows(list(model = model_output, msmt = msmts), .id = "type") %>%
    select(type, time_hr, O2, CO2, DIC) %>%
    pivot_longer(cols = -c(type, time_hr), names_to = "variable") %>%
    left_join(units) %>%
    mutate(label = as.character(label)) %>%
    ggplot() +
    geom_line(aes(x = time_hr / 24,
                  y = value,
                  color = type),
              linewidth = 1.5) +
    facet_wrap(~variable, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("blue", "orange")) +
    # geom_line(data = model_output,
    #           aes(x = time_hr / 24,
    #               y = !!sym(var)),
    #           color = "red") +
    theme_bw() +
    theme(legend.position = c(0.2, 0.8)) +
    labs(x = "time (d)",
         y = expression(mol~m^{-3}))
}


plot_nested_timeseries <- function(df_list, var1, var2) {
  # Merge the dataframes and add a column to identify which dataframe each row comes from
  merged_df <- df_list %>%
    map(~ transform(., dataset = as.character(substitute(.)))) %>% # Get the name of the dataframe as a string
    bind_rows()
  
  # Plot the timeseries with ggplot2
  ggplot(merged_df, aes(x = Date, y = !!sym(var1), color = dataset)) +
    geom_line() +
    geom_line(aes(y = !!sym(var2)), linetype = "dashed") +
    labs(x = "Date", y = "Value", title = "Nested Timeseries Plot")
}
# -------------------------------------
# Author: Jake Diamond
# Purpose: functions to prep measured data for use in the model
# Date: 2023-03-15
# -------------------------------------

# Function to load the data
read_fun <- function(file_path) {
  sheets <- set_names(readxl::excel_sheets(file_path))
  
  # Read each sheet into a data frame
  dfs <- map_dfr(sheets, ~ readxl::read_excel(file_path, sheet = .x) %>%
                   group_by(site) %>% nest(.key = .x) %>% ungroup())
  df <- dfs %>%
    pivot_longer(cols = -site) %>%
    drop_na() %>%
    pivot_wider()
}

# function to get a timeseries from a date range
# also adds a column, "time_hr", which is continuous time in hours
# since the start of the timeseries
ts_f <- function (ts, date_start, date_end) {
  dt <- ts$datetime[2] - ts$datetime[1]
  dt <- as.numeric(dt, units = "hours")
  
  filter(ts, between(date(datetime), date_start, date_end)) %>%
    mutate(time_hr = (row_number() - 1) * dt)
}

# Generates a forcing light signal for the model based on a date range
light_f <- function (data, date_start, date_end, dt = del_t) {
  # Get the x and y values for the approxfun() function
  x <- ts_f(data, date_start, date_end)$time_hr / dt
  y <- ts_f(data, date_start, date_end)$light
  # Create the approxfun() function
  approxfun(x = x,
            y = y,
            method = "linear", rule = 2)
}

# Generates a forcing temp signal for the model based on a date range
temp_f <- function (data, date_start, date_end, dt = del_t) {
  # Get the x and y values for the approxfun() function
  x <- ts_f(data, date_start, date_end)$time_hr / dt
  y <- ts_f(data, date_start, date_end)$temp
  # Create the approxfun() function
  approxfun(x = x, 
            y = y, 
            method = "linear", rule = 2)
}

# Update parameters function
update_parms <- function(x, y, z, replace_vec = "parms") {
  replace_vec <- get(replace_vec)
  dat <- c(x, y, z)
  idx <- intersect(names(dat), names(replace_vec))
  replace_vec[idx] <- dat[idx]
  replace_vec
}

# Update initial conditions functions
update_inits <- function(inis, params, replace_vec = "yini"){
  replace_vec <- get(replace_vec)
  idx <- intersect(names(inis), names(replace_vec))
  replace_vec[idx] <- inis[idx]

  # Equilibrium based on pH and ALK from parameters
  carbeq <- seacarb::carb(flag = 8,
                            var1 = as.numeric(params["pH"]),
                            var2 = as.numeric(params["Alk"]) / 1000,
                            T = as.numeric(params["temp"]),
                            S = 0,
                            warn = "n")

  # Stream DO concentration and carbonate system
  replace_vec$DIC <- carbeq$DIC * 1000        # mol/m3
  replace_vec$ALK <- as.numeric(params["Alk"]) # mol/m3
}

# function to get initial conditions based on date range
init_f <- function (ts, date_start) {
  # First get the inital data from the timeseries if available
  inits <- filter(ts, date(datetime) == date_start) %>%
    head(1)
}

# function to get metabolism/daily conditions based on date range
met_f <- function (data, date_start) {
  filter(data, date == date_start) %>%
    rename(gpp_mean = GPP, er_mean = ER) %>%
    mutate(er_mean = -er_mean)
}

# Function to run the model
mod_fun <- function (parameters, initial_conditions, light, temp) {
  light_signal <<- light
  temp_signal <<- temp
  # Model with user specifications
  ode(y = unlist(initial_conditions),
      times = times,
      func = model,
      parms = unlist(parameters),
      temp_forcing = T,
      light_forcing = T)
}

# function to plot model and measurements
plot_comp_fun <- function (model_output, msmt) {
  p <- ggplot() +
    geom_line(data = msmt,
              aes(x = time_hr / 24,
                  y = O2)) +
    geom_line(data = model_output,
              aes(x = time_hr / 24,
                  y = O2),
              color = "red") +
    theme_bw()
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
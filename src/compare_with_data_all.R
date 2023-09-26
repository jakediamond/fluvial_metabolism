# -------------------------------------
# Author: Jake Diamond
# Purpose: Compare model with measured data
# Date: 2023-04-03
# -------------------------------------
# Load libraries
library(tidyverse)

# Load necessary code -----------------------------------------------------
# source(file.path("src", "carbonate_metabolism_toy_model.R"))
source(file.path("src", "carbonate_model.R"))
source(file.path("src", "compare_functions.R"))
source(file.path("src", "initialize_model.R"))

# Load data ---------------------------------------------------------------
# Filenames of data prepared for comparison
filenames <- list.files(file.path("data"), pattern = "clean.xlsx", 
                        full.names = T, recursive = T)
# Load the data
data <- filenames |>
  map_dfr(read_fun)

# Get the data to model ---------------------------------------------------
# Create a date range for each site
dates <- tribble(
  ~site       , ~datestart  , ~dateend,
  # "FUIROSOS"  , "2019-05-16", "2019-05-17", #small headwater stream in catalonia
  # "REGAS"     , "2019-05-08", "2019-05-09", #small headwater stream in catalonia
  # "SF700"     , "2019-03-07", "2019-03-08", #blackwater river, low ALK in FL
  "SF700"     , "2019-01-21", "2019-01-31", #blackwater river, low ALK in FL
  # "SF2500"    , "2019-06-29", "2019-06-30", #large blackwater river, high ALK in FL
  "SF2500"    , "2019-06-28", "2019-07-10", #large blackwater river, high ALK in FL
  # "ICHE2700"  , "2019-07-01", "2019-07-02", #clearwater spring-fed river in FL
  "ICHE2700"  , "2019-06-03", "2019-06-22", #clearwater spring-fed river in FL
  # "Dampierre" , "2021-08-24", "2021-08-25", #middle Loire River
  # "thom"      , "2019-08-15", "2019-08-16", #Connecticut River
  # "phel"      , "2019-06-25", "2019-06-26", #headwater stream in NE USA
  # "hubb"      , "2019-07-06", "2019-07-07", #Hubbard River in NE USA
  # "nepa"      , "2019-08-29", "2019-08-30", #Nepaug River in NE USA
  # "unio"      , "2019-08-10", "2019-08-11", #Farmington River
) |>
  mutate(
    across(datestart, as.Date),
    across(dateend, as.Date))

# Get the data prepped for the model
mod_data <- data |>
  drop_na() |>
  inner_join(dates) |> # dates to model for each site
  group_by(site) |>
  mutate(
    # time step
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
    mod_inits = map2(inits, mod_parms, update_inits)
  )

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
mods <- mod_data |>
  group_by(site) |>
  transmute(output = pmap(list(mod_parms, mod_inits, dat_comp, delt,
                               lightfun, tempfun), 
                          mod_fun))

# Examine the model output ------------------------------------------------
# Get into good formats and join to data
# mods <- readRDS(file.path("results", "model_comparison_results.RDS"))
mods_out <- transmute(mods, ode_out = map(output, ode_to_df)) |>
  left_join(select(mod_data, dat_comp))
pluck(mods_out, 2, 2)
# Save the data for future use
saveRDS(mods_out, file.path("results", "model_comparison_results_florida.RDS"))

mods_out

plots <- mods_out |>
  group_by(site) |>
  transmute(ps = map2(ode_out, dat_comp, plot_comp_fun))
# plots , var = "O2"
pluck(plots, 2, 3)

patchwork::wrap_plots(plots$ps, ncol = 1, guides = "collect")

pluck(mod_data, "inits", 1)

merged_df <- mods_out |>
  select(-dat_comp) |>
  unnest(ode_out) |>
  pivot_longer(cols = -c(site, time, time_hr))
unnest(dat_comp)
merged_df

ggplot(filter(merged_df, name %in% c("CO2", "O2", "pH", "DIC")),
       aes(x = time_hr / 24,
           y = value)) +
  geom_line() +
  facet_grid(site~name, scales = "free_y")


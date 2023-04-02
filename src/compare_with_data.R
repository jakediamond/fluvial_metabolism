# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
# Load libraries
library(tidyverse)
library(lubridate)

# Load necessary code -----------------------------------------------------
source(file.path("src", "carbonate_metabolism_toy_model.R"))
source(file.path("src", "compare_functions.R"))

# Load data ---------------------------------------------------------------
# Filenames of data prepared for comparison
filenames <- list.files(file.path("data"), pattern = "clean.xlsx", 
                        full.names = T, recursive = T)
# Load the data
data <- filenames %>%
  map_dfr(read_fun)

# Get the data to model ---------------------------------------------------
# Create a date range for each site
dates <- tribble(
  ~site       , ~datestart  , ~dateend,
  "FUIROSOS"  , "2019-05-16", "2019-05-17", #small headwater stream in catalonia
  "REGAS"     , "2019-05-08", "2019-05-09", #small headwater stream in catalonia
  "SF700"     , "2019-03-07", "2019-03-08", #blackwater river, low ALK in FL
  "SF2500"    , "2019-06-29", "2019-06-30", #large blackwater river, high ALK in FL
  "ICHE2700"  , "2019-07-01", "2019-07-02", #clearwater spring-fed river in FL
  "Dampierre" , "2022-05-09", "2022-05-10", #middle Loire River
  "thom"      , "2019-08-16", "2019-08-17", #Connecticut River
  "phel"      , "2019-06-25", "2019-06-26", #headwater stream in NE USA
  "hubb"      , "2019-07-06", "2019-07-07", #Hubbard River in NE USA
  "nepa"      , "2019-08-29", "2019-08-30", #Nepaug River in NE USA
  "unio"      , "2019-08-10", "2019-08-11", #Farmington River
  ) %>%
  mutate(
    across(datestart, as.Date),
    across(dateend, as.Date))

# Get the data prepped for the model
mod_data <- data %>%
  drop_na() %>%
  left_join(dates) %>% # dates to model for each site
  group_by(site) %>%
  mutate(
    # data to compare
    dat_comp = pmap(list(timeseries, datestart, dateend), ts_f),
    # initial conditions of data
    inits = map2(timeseries, datestart, init_f),
    # daily metabolism data
    met = map2(daily, datestart, met_f),
    # light forcing function based on measurements
    lightfun = pmap(list(timeseries, datestart, dateend), light_f),
    # temperature forcing function based on measurements
    tempfun = pmap(list(timeseries, datestart, dateend), temp_f),
    # update model parameters with data
    mod_parms = pmap(list(inits, parms, met), update_parms),
    # update model initial conditions
    mod_inits = map2(inits, mod_parms, update_inits)
  )

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
mods <- mod_data %>%
  # filter(site != "thom") %>%
  group_by(site) %>%
  transmute(output = pmap(list(mod_parms, mod_inits, inits, lightfun, tempfun), 
                          mod_fun))

# Examine the model output ------------------------------------------------
# Get into good formats and join to data
mods_out <- transmute(mods, ode_out = map(output, ode_to_df)) %>%
  left_join(select(mod_data, dat_comp))
mods_out

plots <- mods_out %>%
  group_by(site) %>%
  transmute(ps = map2(ode_out, dat_comp, plot_comp_fun, var = "CO2"))
plots
pluck(plots, 2, 8)

pluck(mod_data, "inits", 1)

mod_data
x=pluck(mod_data, 7, 3)
met_f(x, ymd(20190501))
dplyr::filter(x, date == ymd(20190501)) %>%
  rename(gpp_mean = GPP, er_mean = ER) %>%
  mutate(er_mean = -er_mean)
rm(x)




p
x <- unnest(mods_out, cols = c(ode_out))

pluck(mod_data, 12, 3)

ggplot() +
  geom_line(data = x,
            aes(x = time_hr / 24,
                y = fCO2)) +
  theme_bw() +
  facet_wrap(~site)


y <- pluck(mods_out, 3, 2)
ggplot() +
  geom_line(data = y,
            aes(x = time_hr / 24,
                y = light)) +
  theme_bw()



ggplot(data = x,
       aes(x = time_hr / 24,
           y = LSI)) +
  geom_line() +
  facet_wrap(~site) +
  theme_bw()


x %>%
  mutate(day = ceiling(time_hr / 24)) %>%
  group_by(site, day) %>%
  summarize(gpp = sum(GPP) * 32,
            er = sum(ER) * 32)


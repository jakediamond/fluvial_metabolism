# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------
# Load libraries
library(tidyverse)

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

# Load Loire hourly data
# df_loire <- readRDS(file.path("data", "loire", "all_hourly_data_complete.RDS"))
# Get the data to model ---------------------------------------------------
# Create a date range for each site
dates <- tribble(
  ~site     , ~datestart  , ~dateend,
  "ICHE2700", "2019-07-01", "2019-07-02",
  "FUIROSOS", "2019-05-16", "2019-05-17",
  "REGAS"   , "2019-05-08", "2019-05-09") %>%
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
    mod_parms = pmap(list(parms, inits, met), update_parms),
    # update model initial conditions
    mod_inits = map2(inits, mod_parms, update_inits)
  )

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
mods <- mod_data %>%
  group_by(site) %>%
  transmute(output = pmap(list(mod_parms, mod_inits, lightfun, tempfun), 
                          mod_fun))

# Examine the model output ------------------------------------------------
# Get into good formats and join to data
mods_out <- transmute(mods, ode_out = map(output, ode_to_df)) %>%
  left_join(select(mod_data, dat_comp))
mods_out

p <- mods_out %>%
  mutate(ps = map2(ode_out, dat_comp, plot_comp_fun))

pluck(p, 4, 1)
p
x <- unnest(mods_out, cols = c(ode_out))



ggplot() +
  geom_line(data = x,
            aes(x = time_hr / 24,
                y = O2)) +
  theme_bw() +
  facet_wrap(~site)



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

df_ich %>%
  filter(name %in% c("GPP", "ER")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(day = ceiling(time_hr / 24)) %>%
  group_by(day) %>%
  summarize(gpp = sum(GPP)*32,
            er = sum(ER) * 32)

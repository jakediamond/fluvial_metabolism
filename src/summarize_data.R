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


sum_ts_fun <- function(data) {
  step <- time_length(data$datetime[2] - data$datetime[1], "hour")
  dayl <- 24 / step
  phexists <- "pH" %in% colnames(data)
  dicexists <- "DIC" %in% colnames(data)
  alkexists <- "Alk" %in% colnames(data)
  qexists <- "discharge" %in% colnames(data)
  
  df <- data %>%
    drop_na(O2, CO2) %>%
    summarize(start = min(as.Date(datetime), na.rm = T),
              end = max(as.Date(datetime), na.rm = T),
              O2CO2 = sum(!is.na(O2)) / dayl)
  
  if(phexists & alkexists) {
    co2est <- data %>%
      drop_na(O2, CO2) %>%
      summarize(start = min(as.Date(datetime), na.rm = T),
                end = max(as.Date(datetime), na.rm = T),
                O2CO2est = sum(!is.na(O2)) / dayl)
    df <- co2est
  } else {
    if(phexists) {
      o2co2pH <- data %>%
        drop_na(CO2, O2, pH) %>%
        summarize(O2CO2pH = sum(!is.na(O2)) / dayl) 
      df <- bind_cols(df, o2co2pH)
    }
    if(dicexists) {
      o2co2pHDIC <- data %>%
        drop_na(CO2, O2, pH, DIC) %>%
        summarize(O2CO2pHDIC = sum(!is.na(O2)) / dayl)
      df <- bind_cols(df, o2co2pHDIC)
      
    }
  }
  if(qexists) {
    q <- data %>%
      drop_na(discharge) %>%
      summarize(Q = median(discharge))
    df <- bind_cols(df, q)
  }
  return(df)
}


sum_met_fun <- function(tsdata, metdata) {
  qexists <- "Q" %in% colnames(metdata)
  
  df_dates <- tsdata %>%
    drop_na(O2, CO2) %>%
    mutate(date = date(datetime)) %>%
    distinct(date)
  
  df_met <- semi_join(metdata, df_dates) %>%
    summarize(metcov = n() / nrow(df_dates))
  
  if(qexists) {
    q <- metdata %>%
      drop_na(Q) %>%
      summarize(q = median(Q))
    df_met <- bind_cols(df_met, q)
  }
  
  return(df_met)
}


# pluck(data, 2, 7)
x <- data %>%
  mutate(tssum = map(timeseries, sum_ts_fun),
         metsum = map2(timeseries, daily, sum_met_fun)) %>%
  select(-timeseries, -daily) %>%
  # unnest(metsum)
  unnest(c(tssum, metsum, parms)) %>%
  mutate(Alk = coalesce(Alk, ALK), Q = coalesce(Q, q)) %>%
  select(site, start, end, O2CO2, O2CO2pH, O2CO2pHDIC, O2CO2est, metcov, Q, pH, cond, Ca, Alk, DIC)
write_excel_csv(x, file.path("data", "data_summary.csv"))

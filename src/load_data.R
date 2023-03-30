# -------------------------------------
# Author: Jake Diamond
# Purpose: Load comparison datasets for the model
# Date: 2023-03-15
# -------------------------------------
# Load libraries
library(lubridate)
library(tidyverse)


# Load data ---------------------------------------------------------------
# Load Loire hourly data
df_loire <- readRDS(file.path("data", "loire", "all_hourly_data_complete.RDS"))

# Load Ichetucknee 15-min data
df_ich <- read_csv(file.path("data", "florida", "Ichetucknee", 
                              "ICHE2700_2019_15min_data2.csv")) %>%
  mutate(# make the time local, for some reason it reads in local time
         solar.time = force_tz(solar.time - hours(6)), tz = "US/Eastern")

# Load Santa Fe 15-min data
df_sf <- read_csv(file.path("data", "florida", "Santa Fe", 
                            "SF2500_2019_15min_data.csv")) %>%
  mutate(# make the time local, for some reason it reads in local time
    solar.time = force_tz(solar.time - hours(6)), tz = "US/Eastern")


# read_fun <- function(file_path) {
  sheets <- set_names(readxl::excel_sheets(file_path))
  
  # Read each sheet into a data frame, add a "site" column, and store them in a list
  data_frames <- map(sheets, ~ {
    df <- readxl::read_excel(file_path, sheet = .x)
    # df <- df %>% mutate(site = first(df[["site"]]))
    df
  })
  
  # Combine the data frames into a single data frame, and pivot them wide by sheet
  wide_df <- bind_rows(data_frames, .id = "sheet") %>%
    select(-sheet) %>%
    group_by(site) #%>%
    # pivot_wider(names_from = "sheet", values_from = -c(site))
  
  # Return the wide data frame
  wide_df
# }

test <- read_fun(file.path("data", "florida", "florida_clean.xlsx"))
file_path <- file.path("data", "florida", "florida_clean.xlsx")

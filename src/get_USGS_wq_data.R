# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

# Load the required packages
library(dataRetrieval)
library(tidyverse)

# Define the USGS station codes for the stations of interest
# sites <- c("01187830", "01187300", "01187800", "01188090", "01184000")
sites <- "02320700"
# Set the date range of interest
# start_date <- "-01-01"
# end_date <- "2020-12-31"

# Parameter codes
parameterCd <- c("00400", "00453", "00691", "00915", "00095", "29801", "00688", "39086")
# bicar, ph, dic, 
# Download the water quality data for the stations and date range
data <- readWQPqw(paste0("USGS-", sites), parameterCd)

# Extract the variables of interest and calculate their median values
result <- data %>%
  select(site_no = MonitoringLocationIdentifier, date = ActivityStartDate,
         name = CharacteristicName, value = ResultMeasureValue) %>%
  pivot_wider(names_from = name, values_from = value, values_fn = mean) %>%
  group_by(site_no) %>%
  summarise(across(where(is.numeric), ~median(., na.rm = TRUE)))

# -------------------------------------
# Author: Jake Diamond
# Purpose: Load comparison datasets for the model
# Date: 2023-03-15
# -------------------------------------
# Load libraries
library(tidyverse)


# Function to load the data -----------------------------------------------
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
# Load data ---------------------------------------------------------------
# Filenames of data prepared for comparison
filenames <- list.files(file.path("data"), pattern = "clean.xlsx", 
                        full.names = T, recursive = T)
# Load the data
data <- filenames %>%
  map_dfr(read_fun)

# Load Loire hourly data
df_loire <- readRDS(file.path("data", "loire", "all_hourly_data_complete.RDS"))




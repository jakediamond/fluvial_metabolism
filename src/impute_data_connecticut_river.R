# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

df_thom <- readxl::read_xlsx(file.path("data", "connecticut", 
                                       "connecticut_clean.xlsx"),
                             sheet = "timeseries") %>%
  filter(site == "thom", year(datetime) == 2020)


# Use complete() to fill in missing timesteps and values with NA
df <- df_thom %>% 
  complete(datetime = seq(min(datetime), max(datetime), by = "1 hour"),
           site = "thom")

df_imp <- df %>%
  imputeTS::na_kalman()

ggplot(data = df_imp,
       aes(x = datetime,
           y = DIC)) +
  geom_point()

write_csv(df_imp, file.path("data", "connecticut", "imputed_thom_2020.csv"))

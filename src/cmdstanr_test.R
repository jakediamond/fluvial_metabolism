library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

options(mc.cores = parallel::detectCores())
mod <- cmdstan_model(file.path("src", "carbonate_model.stan"))
?cmdstan_model


parms_data <- unlist(pluck(df_moddat, 13, 1)) # [c("RQ", "PQ", "DIC_l",
#  "ALK_l", "qL")] #"PQ", "O2_l",

parms_data
colnames(df_moddat)
light_signal <- pluck(df_moddat, 11, 1)

temp_signal <- pluck(df_moddat, 12, 1)
ini_data <- unlist(pluck(df_moddat, 14, 1))[c("O2", "DIC", "ALK")]
ini_data
datcomp <- pluck(df_moddat, 8, 1) |>
  select(
    time = time_hr,
    O2,
    CO2,
    pH,
    temp,
    light
  ) |>
  mutate(day = floor(time / 24) + 1) %>%
  as.data.frame()
pHmod <- lm(pH ~ time, data = datcomp)

datcomp$pH <- datcomp$pH - pHmod$coefficients[2] * datcomp$time
plot(datcomp$pH)

yini <- mutate(datcomp, hr = time %% 24) %>%
  filter(hr == 0) %>%
  mutate(c = seacarb::carb(
    flag = 1,
    var1 = pH,
    var2 = CO2 / 1000,
    T = temp,
    S = 0.2,
    pHscale = "F"
  )) %>%
  mutate(
    DIC = c$DIC * 1000,
    ALK = c$ALK * 1000
  ) %>%
  select(O2, DIC, ALK)

data_list <- list(
  N = nrow(datcomp),
  T = 24 / pluck(df_moddat, "delt", 1),
  n_days = as.numeric(difftime(pluck(df_moddat, "dateend", 1),
    pluck(df_moddat, "datestart", 1),
    units = "days"
  )) + 1,
  ts = seq(0, 24 - pluck(df_moddat, "delt", 1), pluck(df_moddat, "delt", 1)),
  y0 = as.matrix(yini),
  y = unlist(split(datcomp[, c("O2", "CO2")], datcomp$day)) |>
    array(dim = c(n_days, 24 / pluck(df_moddat, "delt", 1), 3)),
  forcing =
  )
data_list


fit <- mod$sample(
  data = data_list,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)


fit$summary()
fit$summary(variables = c("theta", "lp__"), "mean", "sd")

# use a formula to summarize arbitrary functions, e.g. Pr(theta <= 0.5)
fit$summary("theta", pr_lt_half = ~ mean(. <= 0.5))

# summarise all variables with default and additional summary measures
fit$summary(
  variables = NULL,
  posterior::default_summary_measures(),
  extra_quantiles = ~ posterior::quantile2(., probs = c(.0275, .975))
)

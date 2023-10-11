
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Simulate data
set.seed(123)
n <- 100
t <- seq(0, 10, length.out = n)
data <- list(n = n, t = t)

# Define ODE model
ode_model <- function(t, y, parms) {
  dydt <- -parms[1] * y[1]
  list(dydt)
}

library(deSolve)
parms <- c(0.1)  # Parameter to be estimated
initial_state <- c(y = 1)  # Initial state
out <- ode(y = initial_state, times = t, func = ode_model, parms = parms)
observed_data <- out[, "y"] + rnorm(n, mean = 0, sd = 0.1)  # Simulated noisy observations

# Fit the model using some method, e.g., least squares
fit <- lm(observed_data ~ t)
plot(observed_data)


library(rstan)

# Compile the Stan model
stan_model <- stan_model(file = file.path("src", 'stan_test_simple.stan'))
# Assuming you have prepared your data in the 'data' list
data_list <- list(
  n = length(observed_data),
  t = data$t,
  y = observed_data
)

# Specify the number of chains and iterations, and other sampling options
chains <- 4
iter <- 2000
warmup <- 1000

# Run the Stan model
fit <- sampling(
  object = stan_model,
  data = data_list,
  chains = chains,
  iter = iter,
  warmup = warmup
)

# Print summary statistics of the posterior
print(fit)

# Display summary statistics
summary(fit)

# Plot trace plots and posterior density plots
traceplot(fit)

stan_trace(fit)
stan_dens(fit, separate_chains = TRUE)
stan_diag(fit)
plot(fit)
# Generate posterior predictive samples
posterior_predictive <- posterior_predict(fit)
get_posterior_mean(fit)
get_
# Plot observed vs. posterior predictive data
plot(data$t, data$y, type = 'l', col = 'blue', ylim = c(min(data$y), max(data$y)))
for (i in 1:chains) {
  lines(data$t, posterior_predictive$data[i, , ], col = 'gray', alpha = 0.5)
}


pred <- as.data.frame(fitNeut, pars = "cPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(xdataNeut)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)
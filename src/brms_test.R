library(brms)
# simulated data for exponential decay: c=A0*exp(-k1*time)
set.seed(42)
time <- c(0,1,2,3,4,5,7,9,11)
y_exp <- 100*exp(-0.2*time)  # c0 =100, kr=0.2, exact simulated data
df_sim <- data.frame(time, y_exp)
df_sim <- df_sim %>% mutate(y_obs=y_exp+rnorm(length(time),0,5)) %>% dplyr::select(-y_exp)# random normal error added to simulated data with mean = 0 and sd = 5
N_t <- length(df_sim$time) # number of datapoints

fo_model2 <- "
  vector fo2(real t,                  //time
             vector y,                // the state
             real k1){                // the parameter next to A0

             vector[1] dydt;          // dimension of ODE
             dydt[1] = - k1 * y[1];   // ODE
             return dydt;             // returns a 1-dimensional vector
  }

// function called from brms: integration of ODEs
vector fo_ode2(data vector time, vector vA0, vector vk1) {
        vector[1] y0 = rep_vector(vA0[1], 1);        //one initial value
        int N = size(time);
        vector[1] y_hat_ode[N] = ode_rk45(fo2, y0 ,0.0,to_array_1d(time),vk1[1]);
        vector[N] y_hat;
        for(i in 1:N) y_hat[i] = y_hat_ode[i,1];
  return(y_hat);
}

"

df_sim_ode <- df_sim[-1,] #remove the first datapoint at t0

fo_formula <- bf(y_obs~fo_ode2(time,A0,k1),
                 A0~1,
                 k1~1,
                 nl=TRUE, loop=FALSE)

fo_priors <- c(prior(normal(100,10),nlpar=A0),
               prior(normal(0.2,0.1), nlpar=k1),
               prior(cauchy(0,10), class=sigma)
)

fo_result2 <- brm(data=df_sim_ode,
                  family = gaussian,
                  formula=fo_formula,
                  prior=fo_priors, 
                  init=0,
                  iter=4000,
                  chains=4, 
                  cores=4, 
                  stanvars = stanvar(scode=fo_model2, 
                                     block="functions"), 
                  backend="cmdstanr", 
                  file="fo_result2")
fo_result2
plot(fo_result2)

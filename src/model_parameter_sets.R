#
# Authors: Jake Diamond
# Purpose: See how model responds to different parameter combos 
# Date: 2023 April 4
# 
# Load functions and model parameters
source(file.path("src", "run_model.R"))
source(file.path("src", "functions_for_updating_model.R"))

# The parameter grid ------------------------------------------------------
# Create dataframe of parameter treatments
trts <- expand.grid(temp = c(10, 15, 20), z = 1, 
                    Alk = c(0.5, 1, 1.5),
                    K600 = 3, #c(1, 2, 5)),
                    gpp_mean = c(10, 12, 15),
                    er_mean = c(10, 12, 15))

# Nice format
df_trts <- mutate(trts, trt = row_number()) %>%
  group_by(trt) %>%
  nest(.key = "trts")

# Get the data prepped for the model
mod_data <- df_trts %>%
  mutate(
    # update model parameters with data
    mod_parms = map(trts, update_parms),
    # update model initial conditions
    mod_inits = map(mod_parms, update_inits)
  )

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
mods <- mod_data %>%
  transmute(output = map2(mod_parms, mod_inits, mod_fun))

# saveRDS(mods, file.path("results", "mod_runs_temp_alk_k600.RDS"))
# saveRDS(mods, file.path("results", "mod_runs_pq_rq_alk_dic_co2.RDS"))
# Look at results ---------------------------------------------------------
# Plot the data
mod_ps <- mods %>%
  mutate(output_df = map(output, ode_to_df),
         plots = map(output_df, plot_fun_all))
trts
x = pluck(mod_ps, 3, 81)
plot(x$DIC)
pluck(mod_ps, 4,29)
trts
plot_fun(pluck(mod_ps, 3, 18), "")
# ggsave(filename = file.path("results", "supp_slope_temp_alk_example_20deg_2alk.png"),
#        dpi = 300,
#        units = "cm",
#        width = 18,
#        height = 12)

trts


p_ex
# ggsave(plot = p_ex,
#        filename = file.path("results", "supp_K_temp_slope_example_20deg_05alk.png"),
#        dpi = 300,
#        units = "cm",
#        width = 18,
#        height = 12)
trts
pluck(mod_ps, 3, 40)
df_comp <- mods %>%
  mutate(output_df = map(output, ode_to_df),
         lms = map(output_df, lm_fun)) %>%
  select(-output, -output_df) %>%
  unnest(cols = lms)
# head(df_comp)

ggplot(data = filter(df_comp, coef == "slope", name == "o2_alk") %>%
         left_join(mutate(trts, trt = row_number())),
       aes(x = temp, y = value, color = Alk)) +
  # scale_color_manual(name = expression(K[600]~"("*d^{-1}*")"), values = c("red", "blue")) +
  geom_point(size = 2) +
  theme_classic(base_size = 10) + 
  scale_color_viridis_c() +
  facet_grid(rows = vars(K600)) +
  # stat_smooth(method = "lm") +
  # ggpubr::stat_regline_equation(size = 3, label.y.npc = "bottom") +
  theme_classic(base_size = 10) +
  # geom_abline(slope = -1) +
  labs(y = expression(O[2]-DIC~slope~"("*beta*")"),
       x = "PQ") + 
  theme(legend.position = c(0.9, 0.5), legend.key.size = unit(0.35, "cm"))


ggplot(data = filter(df_comp, coef == "slope", name == "o2_co2") %>%
         left_join(mutate(trts, trt = row_number())),
       aes(x = PQ, y = value, color = RQ)) +
  # scale_color_manual(name = expression(K[600]~"("*d^{-1}*")"), values = c("red", "blue")) +
  geom_point() +
  theme_classic(base_size = 10) + 
  # facet_grid(rows = vars(PQ), cols = vars(RQ)) +
  # stat_smooth(method = "lm") +
  # ggpubr::stat_regline_equation(size = 3, label.y.npc = "bottom") +
  theme_classic(base_size = 10) +
  geom_abline(slope = -1) +
  labs(y = expression(O[2]-CO[2]~slope~"("*beta*")"),
       x = "PQ") + 
  theme(legend.position = c(0.9, 0.5), legend.key.size = unit(0.35, "cm"))



ggplot(data = filter(df_comp, coef == "slope", name == "o2_co2") %>%
         left_join(mutate(trts, trt = row_number(), pq = gpp_mean/er_mean)),
       aes(x = Alk, y = value, color = as.factor(K600), group = as.factor(K600))) +
  scale_color_manual(name = expression(K[600]~"("*d^{-1}*")"), values = c("red", "blue")) +
  geom_point() +
  theme_classic(base_size = 10) + 
  facet_grid(rows = vars(gpp_mean), cols = vars(er_mean)) +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(size = 3, label.y.npc = "bottom") +
  theme_classic(base_size = 10) +
  labs(y = expression(O[2]-CO[2]~slope~"("*beta*")"),
       x = expression(A[T]~"("*mM*")")) + 
  theme(legend.position = c(0.9, 0.5), legend.key.size = unit(0.35, "cm"))

ggsave(filename = file.path("results", "supp_CO2slope_gpper_alk_k600.png"),
       dpi = 300,
       units = "cm",
       width = 12.2,
       height = 12.2)

p_slope_alk_temp <- ggplot(data = filter(df_comp, coef == "slope", name == "o2_co2") %>%
         left_join(mutate(trts, trt = row_number())),
       aes(x = temp, y = value, color = Alk, group = Alk)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  ggpubr::stat_regline_equation(formula = y ~ x + I(x^2),
                                size = 3, label.y.npc = "bottom") +
  scale_color_viridis_c(name = expression(A[T]~"("*mM*")")) +
  theme_classic(base_size = 10) +
  labs(y = expression(O[2]-CO[2]~slope~"("*beta*")"),
       x = expression(temperature~"("*degree*C*")"))

ggsave(plot = p_slope_alk_temp,
       filename = file.path("results", "supp_slope_temp_alk.png"),
       dpi = 300,
       units = "cm",
       width = 9.2,
       height = 9)




y <- filter(df_comp, coef == "slope", name == "o2_co2") %>%
  rename(temp = trt) %>%
  left_join(x) %>%
  mutate(rat = K_O2 / K_CO2) %>%
  select(temp, slope = value, rat)

p_sloperat <- ggplot(data = y,
       aes(x = rat,
           y = slope,
           color = temp)) +
  geom_point() +
  scale_color_viridis_c(name = expression(temp~"("*degree*C*")")) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.85, 0.4), legend.key.size = unit(0.3, "cm")) + 
  labs(y = expression(O[2]-CO[2]~slope~"("*beta*")"),
       x = expression(K[O[2]]~`:`~K[CO[2]]))




p_slopetemp <- ggplot(data = y,
       aes(x = temp,
           y = slope)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  ggpubr::stat_regline_equation(formula = y ~ x + I(x^2), label.x = 8, label.y = -1.05,
                                size = 3) +
  scale_color_viridis_c() +
  theme_classic(base_size = 10) +
  labs(y = expression(O[2]-CO[2]~slope~"("*beta*")"),
       x = expression(temperature~"("*degree*C*")"))
p_slopetemp


p_all <-p_rat_temp / p_sloperat / p_slopetemp + plot_annotation(tag_levels = "a")
p_all
# ggsave(plot = p_all,
#        filename = file.path("results", "supp_K_temp_slope.png"),
#        dpi = 300,
#        units = "cm",
#        width = 9.2,
#        height = 15)




ggplot(data = filter(df_comp, coef == "slope", name == "o2_alk") %>%
         left_join(mutate(trts, trt = row_number())),
       aes(x = calc, y = value)) +
  # scale_color_manual(name = expression(K[600]~"("*d^{-1}*")"), values = c("red", "blue")) +
  geom_point() +
  theme_classic(base_size = 10) + 
  # facet_grid(rows = vars(gpp_mean), cols = vars(er_mean)) +
  stat_smooth(method = "lm") +
  ggpubr::stat_regline_equation(size = 3, label.y.npc = "bottom") +
  theme_classic(base_size = 10) +
  labs(y = expression(O[2]-CO[2]~slope~"("*beta*")"),
       x = expression(A[T]~"("*mM*")"))

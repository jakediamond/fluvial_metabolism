# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

Sc_CO2 <- Sc("CO2", 1:20)
Sc_O2  <- Sc("O2", 1:20)

# Calculate gas exchange coefficients for CO2 and O2 (1/d)
K_CO2 <- 3 / ((600 / Sc_CO2)^-0.5)
K_O2  <- 3 / ((600 / Sc_O2)^-0.5)
lm(K_CO2~K_O2)
plot(K_CO2, K_O2)

x <- data.frame(K_CO2, K_O2, temp = 1:20)

p_rat_temp <- ggplot(data = x,
       aes(y = K_O2 / K_CO2,
           x = temp)) +
  geom_point() +
  theme_classic(base_size = 10) +
  labs(x = expression(temperature~"("*degree*C*")"),
       y = expression(K[O[2]]~`:`~K[CO[2]]))
p_rat_temp  

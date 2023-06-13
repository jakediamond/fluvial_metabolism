# -------------------------------------
# Author: 
# Purpose: 
# Date:
# -------------------------------------

df <- readRDS(file.path("data", "loire", "all_hourly_data_complete.RDS"))
colnames(df)

df <- df %>%
  mutate(PQ = GPP_mean / ER_mean,
         doy = yday(date),
         regime = if_else(year(date) <2005, "planktonic", "benthic"))

ggplot(data = df,
       aes(x = doy,
           y = PQ,
           color = regime)) +
  stat_smooth()


dftest <- filter(df, date == ymd(20190421))
plot(dftest$DIC)
plot(dftest$O2_mmol_m3)

dftest <- dftest %>%
  mutate(dO2 = O2_mmol_m3 - lag(O2_mmol_m3),
         dDIC = DIC - lag(DIC),
         dALK = (Alk_molkg - lag(Alk_molkg)) * rhow(temp) * 1000,
         TADIC = Alk_molkg * rhow(temp) * 1000 - DIC,
         dCO2 = CO2_uM - lag(CO2_uM),
         ScO2 = Sc("O2", temp),
         ScCO2 = Sc("CO2", temp),
         alpha = sqrt(600/ScO2),
         beta = sqrt(600/ScCO2)) %>%
  mutate(testx = -(O2ex*alpha + exCO2_uM*beta),
         testy = (dO2+dDIC)*24)


ggplot(data = dftest,
       aes(x = testx,
           y = testy,
           color = hr)) +
  geom_point(size = 2) +
  theme_bw() +
  scale_color_viridis_c()+
  ggpubr::stat_regline_equation()



ggplot(data = dftest, #filter(dftest, between(hr, 8, 19)),
       aes(y = O2ex,
           x = CO2ex, #Alk_molkg * rhow(temp) * 1000,
           color = hr)) +
  geom_point(size = 5) +
  theme_bw() +
  scale_color_viridis_c()+
  ggpubr::stat_regline_equation()




df_mod <- df %>%
  select(date, temp, O2 = O2_mmol_m3, CO2 = CO2_uM, ALK = Alk_molkg, DIC) %>%
  mutate(dO2 = O2 - lag(O2),
         dDIC = DIC - lag(DIC),
         dCO2 = CO2 - lag(CO2),
         ALK = ALK * rhow(temp) * 1000,
         dALK = ALK - lag(ALK)) %>%
         # ScO2 = Sc("O2", temp),
         # ScCO2 = Sc("CO2", temp),
         # alpha = sqrt(600/ScO2),
         # beta = sqrt(600/ScCO2)) %>%
  # mutate(testx = -(O2ex*alpha + exCO2_uM*beta),
  #        testy = (dO2+dDIC)*24)
  group_by(date) %>%
  drop_na() %>%
  nest() %>%
  mutate(mod_o2co2 = map(data, ~lm(.$O2~.$CO2)),
         mod_o2dic = map(data, ~lm(.$O2~.$DIC)),
         mod_o2alk = map(data, ~lm(.$O2~.$ALK)),
         mod_alkdic = map(data, ~lm(.$ALK~.$DIC)),
         mod_do2co2 = map(data, ~lm(.$dO2 ~ .$dCO2)),
         mod_do2dic = map(data, ~lm(.$dO2 ~ .$dDIC)))

df_modclean <- df_mod %>%
  mutate(o2co2 = map(mod_o2co2 , broom::tidy),
         o2dic = map(mod_o2dic , broom::tidy),
         o2alk = map(mod_o2alk , broom::tidy),
         alkdic = map(mod_alkdic, broom::tidy),
         do2co2 = map(mod_do2co2, broom::tidy),
         do2dic = map(mod_do2dic, broom::tidy))

df_modsd <- df_modclean %>%
  select(-data, -contains("mod")) %>%
  unnest(c(do2co2, do2dic), names_sep = "_")
  # unnest(c(o2co2, o2dic, o2alk, alkdic), names_sep = "_")

df_mods2d <- df_modsd %>%
  select(date, contains(c("estimate", "term"))) %>%
  pivot_longer(cols = contains("term"), names_sep = "_", names_to = c("model", "term")) %>%
  filter(value != "(Intercept)") %>%
  select(-term, -value) %>%
  pivot_longer(cols = where(is.numeric)) %>%
  select(-model) %>%
  separate_wider_delim(cols = name, delim = "_", names = c("model", "type"))


df_mods2 %>%
  mutate(doy = yday(date)) %>%
  distinct() %>%
  ggplot(aes(x = doy,
             y = value)) +
  # geom_point(alpha = 0.4) +
  stat_summary_bin(bins = 365) +
  theme_bw() +
  facet_wrap(~model, scales = "free_y")

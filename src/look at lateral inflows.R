library(tidyverse)

# Load data ---------------------------------------------------------------
# Filenames of data prepared for comparison
filesq <- list.files(file.path("data", "loire", "discharge"), pattern = ".csv", 
                        full.names = T)
# Load the data
df_q <- vroom::vroom(filesq, .name_repair = janitor::make_clean_names, 
                     id = "filename") %>%
  mutate(site = gsub("^.*?/(.*?)_.*$", "\\1", filename)) %>%
  select(-filename)

# site positions
df_order <- tribble(
  ~site       , ~name           , ~position, ~type,
  "K4080010"  , "loire st satur", 1        , "main",
  "K409401001", "le nohain"     , 2        , "trib",
  "K412301001", "la vrille"     , 3        , "trib",
  "K418001010", "loire gien"    , 4        , "main",
  "K419000201", "la notreure"   , 5        , "trib",
  "K435001020", "loire orleans" , 6        , "main"
)

df_q <- select(df_q, site, date = date_tu, Q = valeur_en_m3_s) %>%
  left_join(df_order)

df_qw <- select(df_q, -site, -name) %>%
  pivot_wider(names_from = c(type, position), values_from = Q)

test <- df_qw %>%
  drop_na(main_1, main_4) %>%
  mutate(qL = main_4 - rowSums(across(c(main_1, trib_2, trib_3)), na.rm = T),
         qL_L = qL / 50115, # distance from st satur to gien (m)
         qL_per = qL / main_4 * 100)
mean(test$qL_L, na.rm = T)
test %>%
  filter(qL_per > 0) %>%
  ggplot(aes(x = month(date),
           y= qL_L)) +
  stat_summary() +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw() +
  labs(x = "month",
       y = "flow fraction from lateral inflow (%)")

ggplot(data = df_qw, 
       aes(x = month(date),
           y = main_4)) +
  stat_summary() +
  scale_x_continuous(breaks = seq(1,12,1)) +
  theme_bw()

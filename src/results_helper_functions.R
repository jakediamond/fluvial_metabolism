# Clean the data
ode_to_df <- function (ode_result) {
  df <- as_tibble(as.data.frame(ode_result))
  colnames(df) <- gsub("\\..*","", colnames(df))
  return(df)
}

# Plotting functions
plot_fun <- function (data, 
                      variables = c("O2", "CO2"), 
                      type = "timeseries",
                      units = units_df) {
  
  # data <- tail(data, 24)
  # Get the units for the x and y variables from the "units_df" dataframe
  var1_lab <- units[units$variable == variables[1], "label"]
  var1_units <- units[units$variable == variables[1], "units"]
  
  if (length(variables) > 1) {
    var2_lab <- units[units$variable == variables[2], "label"]
    var2_units <- units[units$variable == variables[2], "units"]
    if (type == "timeseries" && var1_units == var2_units) {
      var1_lab <- units[units$variable == variables[1], "units_label"]
    }
  }
  
  # If they have the same units and it's a time series, plot together
  
  
  # Determine if plotting CO2 or DIC vs O2 
  abtrue <- ((type == "biplot" && "exCO2" %in% variables && "exO2" %in% variables) ||
               (type == "biplot" && "exDIC" %in% variables && "exO2" %in% variables) ||
               (type == "biplot" && "DIC" %in% variables && "exO2" %in% variables) ||
               (type == "biplot" && "ALK" %in% variables && "exO2" %in% variables) ||
               (type == "biplot" && "ALK" %in% variables && "DIC" %in% variables) ||
             (type == "biplot_daynight" && "exCO2" %in% variables && "exO2" %in% variables) ||
               (type == "biplot_daynight" && "exDIC" %in% variables && "exO2" %in% variables) ||
               (type == "biplot_daynight" && "DIC" %in% variables && "exO2" %in% variables) ||
               (type == "biplot_daynight" && "ALK" %in% variables && "DIC" %in% variables) ||
               (type == "biplot_daynight" && "ALK" %in% variables && "exO2" %in% variables)) 
  
  # Get max and min of datasets
  miny <- min(data[, variables[1]])
  maxy <- max(data[, variables[1]])
  if (length(variables) > 1) {
    minvar2 <- min(data[, variables[2]])
    maxvar2 <- max(data[, variables[2]])
    minx <- miny
    maxx <- maxy
    if (type == "biplot") {
    miny <- minvar2
    maxy <- maxvar2
    } else {
      miny <- min(c(miny, minvar2))
      maxy <- max(c(maxy, maxvar2))
    }
  }
  
  # Do the plotting
  if (type == "timeseries") {
    p <- select(data, time, {{ variables }}) %>%
      pivot_longer(cols = -time) %>%
      ggplot(aes(x = time,
                 y = value,
                 color = name)) +
      geom_line(linewidth = 1.5) +
      theme_bw() +
      geom_hline(yintercept = 0) +
      {if (grepl("sat", variables[1])) list(geom_hline(yintercept = 100))} +
      coord_cartesian(ylim = c(miny, maxy)) +
    labs(x = "time (hr)",
         y = var1_lab)
  }
  if (type == "biplot") {
    var1 <- sym(variables[1])
    var2 <- sym(variables[2])
    p <- ggplot(data = data,
                aes(x = !!var1,
                    y = !!var2,
                    color = time %% 24,
                    group = NULL)) +
      geom_point(size = 2) +
      theme_bw() +
      # geom_hline(yintercept = 0) +
      # geom_vline(xintercept = 0) +
      coord_cartesian(ylim = c(miny, maxy),
                      xlim = c(minx, maxx)) +
      {if (abtrue) list(#geom_abline(slope = -2, intercept = 0, 
                        #            linetype = "dashed"),
                        #geom_abline(slope = -1, intercept = 0, 
                        #            linetype = "dotted"),
                        stat_smooth(method = "lm"),
                        ggpubr::stat_regline_equation(label.x.npc = "center",
                                                      label.y.npc = "center"))} +
      scale_color_viridis_c(name = "hour") +
      labs(x = var1_lab, 
           y = var2_lab)
  }
  if (type == "biplot_daynight") {
    var1 <- sym(variables[1])
    var2 <- sym(variables[2])
    p <- ggplot(data = mutate(data, daynight = if_else(PAR == 0, "night", "day"),
                              day = floor(time / 24)),
                aes(x = !!var1,
                    y = !!var2,
                    color = time %% 24,
                    group = NULL)) +
      geom_point(size = 2) +
      theme_bw() +
      # geom_hline(yintercept = 0) +
      # geom_vline(xintercept = 0) +
      coord_cartesian(ylim = c(miny, maxy),
                      xlim = c(minx, maxx)) +
      facet_grid(rows = vars(day), cols = vars(daynight)) +
      # facet_wrap(~day) +
      {if (abtrue) list(#geom_abline(slope = -2, intercept = 0, 
                        #            linetype = "dashed"),
                        #geom_abline(slope = -1, intercept = 0, 
                        #            linetype = "dotted"),
                        stat_smooth(method = "lm"),
                        ggpubr::stat_regline_equation(label.x.npc = "center",
                                                      label.y.npc = "center"))} +
      scale_color_viridis_c(name = "hour") +
      labs(x = var1_lab, 
           y = var2_lab)
  }
  return(p)
}


# Units dataframe
units_df <- data.frame(
  variable = c("ALK",
               "O2",
               "exO2",
               "O2sat",
               "CO2",
               "exCO2",
               "CO2sat",
               "HCO3",
               "CO3",
               "DIC",
               "exDIC",
               "ER",
               "GPP",
               "calc",
               "fCO2",
               "PAR",
               "SI",
               "pH",
               "temp",
               "GPP_HCO3",
               "LO2",
               "LDIC",
               "LALK"),
  label = I(c(expression(Alkalinity~'('*mol~m^{-3}*')'),
              expression(O[2]~'('*mol~m^{-3}*')'),
              expression(exO[2]~'('*mol~m^{-3}*')'),
              expression(O[2,sat]~'('*`%`~sat.*')'),
              expression(CO[2]~'('*mol~m^{-3}*')'),
              expression(exCO[2]~'('*mol~m^{-3}*')'),
              expression(CO[2,sat]~'('*`%`~sat.*')'),
              expression(HCO[3]^{`-`}~'('*mol~m^{-3}*')'),
              expression(CO[3]^{`2-`}~'('*mol~m^{-3}*')'),
              expression(DIC~'('*mol~m^{-3}*')'),
              expression(exDIC~'('*mol~m^{-3}*')'),
              expression(ER~'('*mol~m^{-2}~h^{-1}*')'),
              expression(GPP~'('*mol~m^{-2}~h^{-1}*')'),
              expression(calcification~'('*mol~m^{-2}~h^{-1}*')'),
              expression(fCO[2]~'('*mol~m^{-2}~h^{-1}*')'),
              expression(PAR~'('*{`mu`}*mol~m^{-2}~s^{-1}*')'),
              expression(SI~'('*`-`*')'),
              expression(pH),
              expression(temperature~'('*degree*C*')'),
              expression(GPP~supported~by~HCO[3]^{`-`}~'('*`%`*')'),
              expression(lateral~O[2]~'('*mol~m^{-2}~h^{-1}*')'),
              expression(lateral~DIC~'('*mol~m^{-2}~h^{-1}*')'),
              expression(lateral~ALK~'('*mol~m^{-2}~h^{-1}*')'))
            ),
  units_label = I(c(expression(mol~m^{-3}),
                    expression(mol~m^{-3}),
                    expression(mol~m^{-3}),
                    expression(`%`~sat.),
                    expression(mol~m^{-3}),
                    expression(mol~m^{-3}),
                    expression(`%`~sat.),
                    expression(mol~m^{-3}),
                    expression(mol~m^{-3}),
                    expression(mol~m^{-3}),
                    expression(mol~m^{-3}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression({`mu`}*mol~m^{-2}~s^{-1}),
                    expression(`-`),
                    expression(`-`),
                    expression(degree*C),
                    expression(`%`),
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}))),
  units = c("mol/m3",
            "mol/m3",
            "mol/m3",
            "%sat",
            "mol/m3",
            "mol/m3",
            "%sat",
            "mol/m3",
            "mol/m3",
            "mol/m3",
            "mol/m3",
            "mol/m2/h",
            "mol/m2/h",
            "mol/m2/h",
            "mol/m2/h",
            "umol/m2/s",
            "-",
            "-",
            "degreeC",
            "%",
            "mol/m2/h",
            "mol/m2/h",
            "mol/m2/h"))

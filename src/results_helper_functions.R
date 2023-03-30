# Clean the data
ode_to_df <- function (ode_result) {
  df <- as_tibble(as.data.frame(ode_result)) %>%
    mutate(time_hr = time * del_t)
  colnames(df) <- gsub("\\..*","", colnames(df))
  return(df)
}

# Plotting functions
plot_fun <- function (data, 
                      variables = c("O2", "CO2"), 
                      type = "timeseries",
                      units = units_df) {
  
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
               (type == "biplot" && "exDIC" %in% variables && "exO2" %in% variables)) 
  
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
    p <- select(data, time_hr, {{ variables }}) %>%
      pivot_longer(cols = -time_hr) %>%
      ggplot(aes(x = time_hr,
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
                    color = time_hr %% 24,
                    group = NULL)) +
      geom_point(size = 2) +
      theme_bw() +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      coord_cartesian(ylim = c(miny, maxy),
                      xlim = c(minx, maxx)) +
      {if (abtrue) list(geom_abline(slope = -2, intercept = 0, 
                                    linetype = "dashed"),
                        geom_abline(slope = -1, intercept = 0, 
                                    linetype = "dotted"),
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
               "DIC",
               "exDIC",
               "ER",
               "GPP",
               "calc",
               "fCO2",
               "PAR",
               "LSI",
               "pH",
               "temp"),
  label = I(c(expression(Alkalinity~'('*mol~m^{-3}*')'),
              expression(O[2]~'('*mol~m^{-3}*')'),
              expression(exO[2]~'('*mol~m^{-3}*')'),
              expression(O[2,sat]~'('*`%`~sat.*')'),
              expression(CO[2]~'('*mol~m^{-3}*')'),
              expression(exCO[2]~'('*mol~m^{-3}*')'),
              expression(CO[2,sat]~'('*`%`~sat.*')'),
              expression(DIC~'('*mol~m^{-3}*')'),
              expression(exDIC~'('*mol~m^{-3}*')'),
              expression(ER~'('*mol~m^{-2}~h^{-1}*')'),
              expression(GPP~'('*mol~m^{-2}~h^{-1}*')'),
              expression(calcification~'('*mol~m^{-2}~h^{-1}*')'),
              expression(fCO[2]~'('*mol~m^{-2}~h^{-1}*')'),
              expression(PAR~'('*{`mu`}*mol~m^{-2}~s^{-1}*')'),
              expression(LSI~'('*`-`*')'),
              expression(pH),
              expression(temperature~'('*degree*C*')'))
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
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression(mol~m^{-2}~h^{-1}),
                    expression({`mu`}*mol~m^{-2}~s^{-1}),
                    expression(`-`),
                    expression(`-`),
                    expression(degree*C))),
  units = c("mol/m3",
            "mol/m3",
            "mol/m3",
            "%sat",
            "mol/m3",
            "mol/m3",
            "%sat",
            "mol/m3",
            "mol/m3",
            "mol/m2/h",
            "mol/m2/h",
            "mol/m2/h",
            "mol/m2/h",
            "umol/m2/s",
            "-",
            "-",
            "degreeC"))

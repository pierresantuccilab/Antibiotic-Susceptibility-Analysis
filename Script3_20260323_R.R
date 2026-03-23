############################################################
# SCRIPT NAME: DoseResponse_TwoStrains_Comparison.R
#
# PURPOSE
# This script performs dose-response curve (DRC) analysis
# for two bacterial strains exposed to a given compound.
#
# It:
#   1. Filters the dataset for selected compound and strains
#   2. Cleans and formats the data
#   3. Computes mean growth per replicate
#   4. Fits dose-response models (4-parameter logistic)
#   5. Plots the curves
#   6. Calculates EC50 and EC90 values
#
# INPUT DATA REQUIREMENTS
# Dataframe: Dataset obtained after Step2 and save as CSV file
# Required columns:
#   - compound         : compound name (character)
#   - strain           : strain name (character)
#   - concentrations   : compound concentration (numeric or character)
#   - perodc           : relative growth (%) (numeric or character)
#   - repb             : biological replicate ID
#
# OUTPUT
# - Dose-response plot (R plot window / exportable)
# - EC50 and EC90 values (displayed in plot legend)
#
# ASSUMPTIONS
# - perodc values are normalized (0–100 or higher)
# - Negative values are artifacts and set to 0
# - Two strains are compared simultaneously
#
# USER-DEFINED PARAMETERS (EDIT THESE)
compound_name <- "Genta"
strain_1 <- "WT"
strain_2 <- "AmpR"

x_label <- "Genta (mg/L)"
plot_title <- "Dose-response: WT vs AmpR"

x_limits <- c(0, 64)
y_limits <- c(0, 200)

color_1 <- "#67A9CF"
color_2 <- "#FDD0A2"
############################################################


############################
# STEP 1 — Load packages
############################
library(tidyverse)
library(dplyr)
library(drc)

options(scipen = 999)  # Disable scientific notation


############################
# STEP 2 — Filter dataset
############################
data <- dplyr::filter(
  Data_Compilation_Strain1vsStrain2,
  compound == "Genta",
  strain %in% c(strain_1, strain_2)
)


############################
# STEP 3 — Data cleaning
############################
data$concentrations <- as.numeric(as.character(data$concentrations))
data$perodc <- as.numeric(as.character(data$perodc))

# Replace negative values with 0
data$perodc[data$perodc < 0] <- 0


############################
# STEP 4 — Compute mean per replicate
############################
data_mean <- data %>%
  group_by(strain, concentrations, repb) %>%
  summarise(
    perodc_mean = mean(perodc, na.rm = TRUE),
    .groups = "drop"
  )


############################
# STEP 5 — Fit dose-response models
############################
fit_1 <- drm(
  perodc_mean ~ concentrations,
  data = filter(data_mean, strain == strain_1),
  fct = LL.4()
)

fit_2 <- drm(
  perodc_mean ~ concentrations,
  data = filter(data_mean, strain == strain_2),
  fct = LL.4()
)


############################
# STEP 6 — Plot curves
############################
plot(fit_1,
     broken = TRUE,
     type = "all",
     xlab = x_label,
     ylab = "Relative growth (%)",
     main = plot_title,
     xlim = x_limits,
     ylim = y_limits,
     lwd = 2,
     pch = 21, cex=2,
     bg = color_1,
     col = color_1)

plot(fit_2,
     broken = TRUE,
     type = "all",
     add = TRUE,
     lwd = 2,
     pch = 21, cex=2,
     bg = color_2,
     col = color_2)


############################
# STEP 7 — Calculate EC values
############################
ec50_1 <- ED(fit_1, 50, interval = "delta")
ec50_2 <- ED(fit_2, 50, interval = "delta")

EC50_1_val <- ec50_1[1, 1]
EC50_1_se  <- ec50_1[1, "Std. Error"]

EC50_2_val <- ec50_2[1, 1]
EC50_2_se  <- ec50_2[1, "Std. Error"]


ec90_1 <- ED(fit_1, 90, interval = "delta")
ec90_2 <- ED(fit_2, 90, interval = "delta")

EC90_1_val <- ec90_1[1, 1]
EC90_1_se  <- ec90_1[1, "Std. Error"]

EC90_2_val <- ec90_2[1, 1]
EC90_2_se  <- ec90_2[1, "Std. Error"]


############################
# STEP 8 — Add legends
############################
legend("topright",
       legend = c(strain_1, strain_2),
       col = c(color_1, color_2),
       pch = 21,
       pt.bg = c(color_1, color_2),
       lwd = 2)

legend("topleft",
       legend = c(
         bquote(EC[50]~"(mg/L)"),
         bquote(.(strain_1)~":"~.(round(EC50_1_val, 2))~"\u00B1"~.(round(EC50_1_se, 2))),
         bquote(.(strain_2)~":"~.(round(EC50_2_val, 2))~"\u00B1"~.(round(EC50_2_se, 2)))
       ),
       cex = 0.8)

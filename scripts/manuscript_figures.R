#--------------------------------------------------------------------------- ---
# Code to produce figures for Bratt et al
# Client: Proteus
# Year: 2024
#--------------------------------------------------------------------------- ---

#--------------------------------------------------------------------------- ---
# Version update history:
#   Abby Bratt (18/11/24)
#--------------------------------------------------------------------------- ---

## SET WORKING DIRECTORIES ----

# Point to code repository
repo_dir <- getwd()

#Point to data repository
SP_dir <- ProteusFunctions::find_OD_dir("PPM")

## LOAD LIBRARIES ----

library(tidyverse)
library(lubridate)
library(bayesplot)
library(nimble)
library(here)
library(beepr)
library(postpack)
library(pbapply)
library(strex)
library(RColorBrewer)
library(cowplot)
library(ProteusFunctions) # company functions
library(ProteusTheme) # company ggplot theme

# SHORT TIMESCALE ----

## Load data ----

# Load processed density and occupancy data
load(here("Data", "short_pre_jags_image.RData"))

# Load MCMC results
load(here("Results", "sqrt", "reanalysis_short_samples_sqrt.Rdata"))
load(here("Results", "sqrt", "reanalysis_short_summary_sqrt.Rdata"))
est <- read_csv(here("Results", "sqrt", "short time density and occupancy estimates_reanalysis_sqrt.csv")) %>%
  mutate(`100mPlot` = as_factor(`100mPlot`))

# Load predictions
AbundanceEstimates_3 <- read_csv(here("Results", "sqrt", "reanalysis_predictions_short_sqrt_upscaled.csv")) %>%
  mutate(my = as.Date(paste(Year, as.numeric(Month_fct) + 3, "01", sep = "-"))) %>%
  filter(Year < 2022) %>%
  mutate(
    Site = str_to_sentence(Site),
    Site = if_else(Site == "Ssm", "SSM", Site)
  )

## Density vs occupancy ----

est <- est %>%
  mutate(
    Site = str_to_sentence(Site),
    Site = if_else(Site == "Ssm", "SSM", Site)
  )

sp1 <- ggplot(est, aes(x = prop_occ_mean, y = den_mean, color = `Site`)) +
  geom_pointrange(aes(xmin = prop_occ_l_90, xmax = prop_occ_u_90), colour = "grey", alpha = 0.35, size = 0.25) +
  geom_pointrange(aes(ymin = den_l_90, ymax = den_u_90), colour = "grey", alpha = 0.35, size = 0.25) +
  geom_point(shape = 19, size = 1.5, alpha = 0.9) +
  ylab("Estimated density") +
  xlab("Estimated PSO") +
  theme_proteus() +
  theme(panel.grid.major.y = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_color_brewer(palette = "Dark2")
sp1

ggsave(here("Results", "ms", "MS_reanalysis_short_est_den_occ_sqrt.png"),
  width = 15, height = 15, units = "cm", dpi = 320
)

## Density vs detection ----

sp2 <- ggplot(est, aes(x = occ_p_mean, y = den_mean, color = `Site`)) +
  geom_pointrange(aes(xmin = occ_p_l_90, xmax = occ_p_u_90), colour = "grey", alpha = 0.35, size = 0.25) +
  geom_pointrange(aes(ymin = den_l_90, ymax = den_u_90), colour = "grey", alpha = 0.35, size = 0.25) +
  geom_point(shape = 19, size = 1.5, alpha = 0.9) +
  ylab("Estimated density") +
  xlab("Estimated detection probability") +
  theme_proteus() +
  theme(panel.grid.major.y = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_color_brewer(palette = "Dark2")
sp2

ggsave(here("Results", "ms", "MS_reanalysis_short_est_den_occ_p_sqrt.png"),
  width = 15, height = 15, units = "cm", dpi = 320
)

## Density timeseries ----

est <- est %>%
  mutate(my = as.Date(paste(Year, as.numeric(Month_fct) + 3, "01", sep = "-")))

sp3 <- ggplot(est, aes(x = my, y = den_mean, color = `Site`)) +
  geom_errorbar(aes(ymin = den_l_90, ymax = den_u_90),
    size = 0.25, width = 0, position = position_dodge(width = 1)
  ) +
  geom_point(shape = 19, size = 1, alpha = 0.9, position = position_dodge(width = 1)) +
  facet_wrap(~Site, ncol = 3) +
  ylab("Estimated density") +
  theme_proteus() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Year") +
  theme(strip.background = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0))
sp3

ggsave(here("Results", "ms", "MS_reanalysis_short_den_year_sqrt.png"),
  width = 20, height = 10, units = "cm", dpi = 320
)

## Abundance timeseries ----

sp4 <- ggplot(AbundanceEstimates_3, aes(x = my, y = mean_N, colour = Site)) +
  geom_errorbar(aes(ymin = l_90_N, ymax = u_90_N),
    size = 0.25, width = 0, position = position_dodge(width = 1)
  ) +
  geom_point(shape = 19, size = 1, alpha = 0.9, position = position_dodge(width = 1)) +
  facet_wrap(~Site, ncol = 3) +
  ylab("Estimated abundance") +
  theme_proteus() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Year") +
  theme(strip.background = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0))
sp4

ggsave(here("Results", "ms", "MS Scaled abundance estimates SHORT.png"),
  width = 20, height = 15, units = "cm", dpi = 320
)

# LONG TIMESCALE ----

## Load data ----

# Load processed density and occupancy data
load(here("Data", "long_prejags_image.RData"))

# Load MCMC results
load(here("Results", "sqrt", "reanalysis_long_samples_sqrt.Rdata"))
load(here("Results", "sqrt", "reanalysis_long_summary_sqrt.Rdata"))
est <- read_csv(here("Results", "sqrt", "long time density and occupancy estimates_reanalysis_sqrt.csv")) %>%
  mutate(`100mPlot` = as_factor(`100mPlot`))

# Load predictions
AbundanceEstimates_3 <- read_csv(here("Results", "sqrt", "reanalysis_predictions_sqrt_upscaled.csv")) %>%
  filter(Year < 2022) %>%
  mutate(
    Site = str_to_sentence(Site),
    Site = if_else(Site == "Ssm", "SSM", Site)
  )

## Density vs occupancy ----

lp1 <- ggplot(est, aes(x = prop_occ_mean, y = den_mean, color = `Site`)) +
  geom_pointrange(aes(xmin = prop_occ_l_90, xmax = prop_occ_u_90),
    colour = "grey", alpha = 0.35, size = 0.25
  ) +
  geom_pointrange(aes(ymin = den_l_90, ymax = den_u_90),
    colour = "grey", alpha = 0.35, size = 0.25
  ) +
  geom_point(shape = 19, size = 1.5, alpha = 0.9) +
  ylab("Estimated density") +
  xlab("Estimated PSO") +
  theme_proteus() +
  theme(panel.grid.major.y = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_color_brewer(palette = "Dark2")
lp1

ggsave(here("Results", "ms", "MS_reanalysis_long_est_den_occ_sqrt.png"),
  width = 15, height = 15, units = "cm", dpi = 320
)

## Density vs detection ----

lp2 <- ggplot(est, aes(x = occ_p_mean, y = den_mean, color = `Site`)) +
  geom_pointrange(aes(xmin = occ_p_l_90, xmax = occ_p_u_90),
    colour = "grey", alpha = 0.35, size = 0.25
  ) +
  geom_pointrange(aes(ymin = den_l_90, ymax = den_u_90),
    colour = "grey", alpha = 0.35, size = 0.25
  ) +
  geom_point(shape = 19, size = 1.5, alpha = 0.9) +
  ylab("Estimated density") +
  xlab("Estimated detection probability") +
  theme_proteus() +
  theme(panel.grid.major.y = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_color_brewer(palette = "Dark2")
lp2

ggsave(here("Results", "ms", "MS_reanalysis_long_est_den_occ_p_sqrt.png"),
  width = 15, height = 15, units = "cm", dpi = 320
)

## Density timeseries ----

lp3 <- ggplot(est, aes(x = Year, y = den_mean, color = `Site`)) +
  geom_errorbar(aes(ymin = den_l_90, ymax = den_u_90),
    size = 0.25, width = 0, position = position_dodge(width = 1)
  ) +
  geom_point(shape = 19, size = 1, alpha = 0.9, position = position_dodge(width = 1)) +
  facet_wrap(~Site, ncol = 3) +
  ylab("Estimated density") +
  theme_proteus() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Year") +
  theme(strip.background = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0))
lp3

ggsave(here("Results", "ms", "MS_reanalysis_long_den_year_sqrt.png"),
  width = 20, height = 10, units = "cm", dpi = 320
)

## Abundance timeseries ----

lp4 <- ggplot(AbundanceEstimates_3, aes(x = Year, y = mean_N, colour = Site)) +
  geom_errorbar(aes(ymin = l_90_N, ymax = u_90_N),
    size = 0.25, width = 0, position = position_dodge(width = 1)
  ) +
  geom_point(shape = 19, size = 1, alpha = 0.9, position = position_dodge(width = 1)) +
  facet_wrap(~Site, ncol = 3) +
  ylab("Estimated abundance") +
  theme_proteus() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Year") +
  theme(strip.background = element_blank()) +
  scale_y_continuous(expand = c(0.01, 0))
lp4

ggsave(here("Results", "ms", "MS Scaled abundance estimates LONG.png"),
  width = 20, height = 15, units = "cm", dpi = 320
)


# COMBINED PLOTS ----

plot_grid(sp1, lp1, labels = c("A", "B"), nrow = 1, ncol = 2)
ggsave(here("Results", "ms", "MS_reanalysis_both_est_den_occ_sqrt.png"),
  width = 30, height = 15, units = "cm", dpi = 320
)

plot_grid(sp2, lp2, labels = c("A", "B"), nrow = 1, ncol = 2)
ggsave(here("Results", "ms", "MS_reanalysis_both_est_den_occ_p_sqrt.png"),
  width = 30, height = 15, units = "cm", dpi = 320
)

plot_grid(sp3, lp3, labels = c("A", "B"), nrow = 2, ncol = 1)
ggsave(here("Results", "ms", "MS_reanalysis_both_den_year_sqrt.png"),
  width = 20, height = 20, units = "cm", dpi = 320
)

foo <- plot_grid(sp4, lp4, labels = c("A", "B"), nrow = 2, ncol = 1)

ggdraw() +
  draw_plot(foo) + # Overlay PPM artwork onto plot
  draw_image(
    image = here("Data", "20220328_ppm_color_transparent.png"),
    x = 0.925, y = 0.85, hjust = 0.5, vjust = 0.5,
    halign = 0.5, valign = 0.5, scale = 0.125
  ) +
  draw_image(
    image = here("Data", "20220328_ppm_color_transparent.png"),
    x = 0.925, y = 0.35, hjust = 0.5, vjust = 0.5,
    halign = 0.5, valign = 0.5, scale = 0.125
  )

ggsave(here("Results", "ms", "MS Scaled abundance estimates BOTH.png"),
  width = 20, height = 20, units = "cm", dpi = 320
)

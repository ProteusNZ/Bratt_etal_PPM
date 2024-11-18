#--------------------------------------------------------------------------- ---
# Code to produce in-text results for Bratt et al
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

# Point to data repository
SP_dir <- ProteusFunctions::find_OD_dir("USGS")

## LOAD LIBRARIES ----

library(here)
library(tidyverse)
library(postpack)
library(nimble)

## LOAD DATA ----

### Short timescale ----

# Model fit
load(here("Results/sqrt/reanalysis_short_summary_sqrt.Rdata"))
summ.short <- summ

# MCMC samples
load(here("Results/sqrt/reanalysis_short_samples_sqrt.Rdata"))
jags.short <- jagsfit

# Predictions 
est.short <-read_csv(here("Results","sqrt", 
                          "short time density and occupancy estimates_reanalysis_sqrt.csv")) %>% 
  mutate(`100mPlot` = as_factor(`100mPlot`))
abund.short <-read_csv(here("Results","sqrt", 
                            "reanalysis_predictions_short_sqrt_upscaled.csv")) %>% 
  filter(Year < 2022)

### Long timescale ----

# Model fit
load(here("Results/sqrt/reanalysis_long_summary_sqrt.Rdata"))
summ.long <- summ

# MCMC samples
load(here("Results/sqrt/reanalysis_long_samples_sqrt.Rdata"))
jags.long <- jagsfit

# Predictions 
est.long <-read_csv(here("Results","sqrt", 
                         "long time density and occupancy estimates_reanalysis_sqrt.csv")) %>% 
  mutate(`100mPlot` = as_factor(`100mPlot`))
abund.long <-read_csv(here("Results","sqrt", 
                           "reanalysis_predictions_sqrt_upscaled.csv")) %>% 
  filter(Year < 2022)

rm(summ, jagsfit)

# DENSITY ----

## Short timescale ----

# expected density on core plots
est.short %>% 
  summarise(mean = mean(den_mean, na.rm = T))
est.short %>% 
  summarise(mean = quantile(den_mean, 0.025, na.rm = T))
est.short %>% 
  summarise(mean = quantile(den_mean, 0.975, na.rm = T))

# effective area sampled
pars <- post_subset(jags.short, c("A"))
post_summ(pars, "A")

# mean capture probability - single check
pars <- post_subset(jags.short, c("mu_hug_p"))
post_summ(pars, "mu_hug_p") %>% 
  plogis()

# cumulative capture probability - 6 checks
1 - (1 - (post_summ(pars, "mu_hug_p") %>% plogis()))^6

## Long timescale ----

# expected density on core plots
est.long %>% 
  summarise(mean = mean(den_mean, na.rm = T))
est.long %>% 
  summarise(mean = quantile(den_mean, 0.025, na.rm = T))
est.long %>% 
  summarise(mean = quantile(den_mean, 0.975, na.rm = T))

# effective area sampled
pars <- post_subset(jags.long, c("A"))
post_summ(pars, "A")

# mean capture probability - single check
pars <- post_subset(jags.long, c("mu_hug_p"))
post_summ(pars, "mu_hug_p") %>% 
  plogis()

# cumulative capture probability - 4 checks
1 - (1 - (post_summ(pars, "mu_hug_p") %>% plogis()))^4

# OCCUPANCY ----

## Short timescale ----

# expected occupancy on all plots
est.short %>% 
  summarise(mean = mean(psi_mean, na.rm = T))
est.short %>% 
  summarise(mean = quantile(psi_mean, 0.025, na.rm = T))
est.short %>% 
  summarise(mean = quantile(psi_mean, 0.975, na.rm = T))

# detection
est.short %>% 
  summarise(mean = mean(occ_p_mean, na.rm = T))
est.short %>% 
  summarise(mean = quantile(occ_p_mean, 0.025, na.rm = T))
est.short %>% 
  summarise(mean = quantile(occ_p_mean, 0.975, na.rm = T))

# cumulative capture probability - 2 checks
1 - (1 - (
  est.short %>% 
            summarise(mean = mean(occ_p_mean, na.rm = T))
  ))^2
1 - (1 - (
  est.short %>% 
    summarise(mean = quantile(occ_p_mean, 0.025, na.rm = T))
))^2
1 - (1 - (
  est.short %>% 
    summarise(mean = quantile(occ_p_mean, 0.975, na.rm = T))
))^2

# psi grid type effect 
# detection grid type effect
cov.pars <- post_subset(jags.short, c("beta_"))
post_summ(cov.pars, get_params(cov.pars))

## Long timescale ----

# expected occupancy on all plots
est.long %>% 
  summarise(mean = mean(psi_mean, na.rm = T))
est.long %>% 
  summarise(mean = quantile(psi_mean, 0.025, na.rm = T))
est.long %>% 
  summarise(mean = quantile(psi_mean, 0.975, na.rm = T))

# detection
est.long %>% 
  summarise(mean = mean(occ_p_mean, na.rm = T))
est.long %>% 
  summarise(mean = quantile(occ_p_mean, 0.025, na.rm = T))
est.long %>% 
  summarise(mean = quantile(occ_p_mean, 0.975, na.rm = T))

# cumulative capture probability - 2 checks
1 - (1 - (
  est.long %>% 
    summarise(mean = mean(occ_p_mean, na.rm = T))
))^7
1 - (1 - (
  est.long %>% 
    summarise(mean = quantile(occ_p_mean, 0.025, na.rm = T))
))^7
1 - (1 - (
  est.long %>% 
    summarise(mean = quantile(occ_p_mean, 0.975, na.rm = T))
))^7

# psi grid type effect 
# detection grid type effect
cov.pars <- post_subset(jags.long, c("beta_"))
post_summ(cov.pars, get_params(cov.pars))

# PREDICTED ABUNDANCE ----

## Short timescale ----

plot_dat <- abund.short %>% 
  select(Site, Year, Month_fct, mean_N, sd_N,
         contains("90")) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 0))) %>% 
  mutate(crI90 = paste0("(", `l_90_N`, ", ", `u_90_N`, ")")) %>% 
  mutate(Month = month.abb[Month_fct+3], .before = 4) %>% 
  select(-c(3, 7:8))

# across all sites, predicted N
plot_dat %>% 
  group_by(Site, Month) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  select(-Year) %>% 
  arrange(desc(mean_N))

# across all years, predicted N
plot_dat %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  select(-Year) %>% 
  arrange(desc(mean_N))

# across all years, density
est.short %>% 
  select(Site, `100mPlot`, Year, month,
         contains("den")&!contains("5")&!contains("med")) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 2)), 
         month = month.abb[month]) %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  select(-Year) %>% 
  arrange(desc(den_mean))

# average population growth rate - all sites
plot_dat <- abund.short %>% 
  select(Site, Year, Month_fct, mean_N, sd_N,
         contains("90")) %>% 
  mutate(Month = month.abb[Month_fct+3], .before = 4) %>% 
  select(-Month_fct) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 0))) %>% 
  group_by(Site, Year) %>% 
  mutate(lambda = mean_N/lag(mean_N)) %>% 
  ungroup() %>% 
  summarise(mean_lambda  = mean(lambda, na.rm = T))
plot_dat

# average population growth rate - by site
plot_dat <- abund.short %>% 
  select(Site, Year, Month_fct, mean_N, sd_N,
         contains("90")) %>% 
  mutate(Month = month.abb[Month_fct+3], .before = 4) %>% 
  select(-Month_fct) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 0))) %>% 
  group_by(Site, Year) %>% 
  mutate(lambda = mean_N/lag(mean_N)) %>% 
  summarise(mean_lambda  = mean(lambda, na.rm = T))
plot_dat

## Long timescale ----

plot_dat <- abund.long %>% 
  select(Site, Year, mean_N, sd_N,
         contains("90")) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 0))) 

# across all sites, predicted N
plot_dat %>% 
  group_by(Year) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  arrange(desc(mean_N))

# across all years, predicted N
plot_dat %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  select(-Year) %>% 
  arrange(desc(mean_N))

# across all years, density
est.long %>% 
  select(Site, `100mPlot`, Year, 
         contains("den")&!contains("5")&!contains("med")) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 2))) %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  select(-Year) %>% 
  arrange(desc(den_mean))

# average population growth rate - all sites
plot_dat <- abund.long %>% 
  select(Site, Year, mean_N, sd_N,
         contains("90")) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 0))) %>% 
  group_by(Site) %>% 
  mutate(lambda = mean_N/lag(mean_N)) %>% 
  ungroup() %>% 
  summarise(mean_lambda  = mean(lambda, na.rm = T))
plot_dat

# average population growth rate - by site
plot_dat <- abund.long %>% 
  select(Site, Year, mean_N, sd_N,
         contains("90")) %>% 
  mutate(across(where(is.numeric), ~round(., digits = 0))) %>% 
  group_by(Site) %>% 
  mutate(lambda = mean_N/lag(mean_N)) %>% 
  summarise(mean_lambda  = mean(lambda, na.rm = T))
plot_dat


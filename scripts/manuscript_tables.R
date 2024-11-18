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
library(rjags)
library(here)
library(beepr)
library(postpack)
library(pbapply) 
library(strex)
library(RColorBrewer)
library(ProteusFunctions)
library(ProteusTheme)

# SHORT TIMESCALE ----

## Load data ----

# Load processed density and occupancy data
load(here("Data", "short_pre_jags_image.RData"))

# Load MCMC results
load(here("Results", "sqrt", "reanalysis_short_samples_sqrt.Rdata"))
load(here("Results", "sqrt", "reanalysis_short_summary_sqrt.Rdata"))

# Load predictions
est <-read_csv(here("Results","sqrt", "short time density and occupancy estimates_reanalysis_sqrt.csv")) %>% 
  mutate(`100mPlot` = as_factor(`100mPlot`))

## Create tables ----

model.pars <- do.call(rbind, model.pars) %>% as.data.frame() %>% 
  pivot_longer(everything())

model_sel <- model.pars %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), 
            sd = sd(value), 
            l_95 = quantile(value, 0.025), 
            l_90 = quantile(value, 0.05), 
            l_50 = quantile(value, 0.25), 
            med = quantile(value, 0.5),
            u_50 = quantile(value, 0.75),
            u_90 = quantile(value, 0.95),
            u_95 = quantile(value, 0.975)
  ) %>% 
  mutate(index = case_when(
    is.na(str_first_number(name)) ~ "Random",
    str_first_number(name) == 1 ~ "Intercept",
    str_first_number(name) == 2 ~ "PSO effect",
    str_first_number(name) == 3 ~ "Detection effect"
  ),
  index=factor(index,levels=c("Intercept","PSO effect","Detection effect","Random")),
  .after = 1) %>% 
  mutate(name = case_when(
    str_detect(name, "beta") & index == "" ~ "Base", 
    str_detect(name, "beta") ~ "Base", 
    str_detect(name, "sd") ~ "SD", 
  )) %>% 
  rename(label = name)



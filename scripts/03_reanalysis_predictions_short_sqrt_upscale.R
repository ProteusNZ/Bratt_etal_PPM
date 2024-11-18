#--------------------------------------------------------------------------- ---
# Code to upscale abundance, short timescale
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

## LOAD DATA SETS ----

### Preprocessed data ----
load(here("Data", "short_pre_jags_image.RData"))
load(here("Data", "corrected_X.rda"))

SiYr_wData <- SiYr_all[sel, ]
SiYr_wData <- SiYr_wData %>%
  mutate(Year = as.numeric(as.character(Year))) %>%
  left_join(X %>% select(Year, Site, `100mPlot`, SiteType, SiteTypeCode) %>% distinct(),
    by = c("Site", "Year", "100mPlot")
  )

### Model results ----

load(here("reanalysis_long_samples_sqrt.Rdata"))

A.estimates <- do.call(rbind, post_subset(jagsfit, "A")) %>%
  as.data.frame() %>%
  mutate(index = row_number())

N.estimates2 <- do.call(rbind, post_subset(jagsfit, "^N")) %>%
  as.data.frame() %>%
  pivot_longer(everything()) %>%
  mutate(index = str_first_number(name)) %>%
  left_join(SiYr_wData %>% rownames_to_column(var = "index") %>% mutate(index = as.numeric(index)),
    by = "index"
  ) %>%
  # filter(SiteType != "PermCore") %>% #####
  mutate(TotalSiteArea = case_when(
    (Site == "EDSON" & SiteTypeCode == 1) ~ 471,
    (Site == "EDSON" & SiteTypeCode == 2) ~ 3,
    (Site == "OSCAR" & SiteTypeCode == 1) ~ 407,
    (Site == "OSCAR" & SiteTypeCode == 2) ~ 4,
    (Site == "SSM" & SiteTypeCode == 1) ~ 99,
    (Site == "SSM" & SiteTypeCode == 2) ~ 6
  ))

## ESTIMATE ABUNDANCE ----

iter_idx <- sample(1:90000, 10000)

AbundanceEstimates_3 <- N.estimates2 %>%
  group_by(Site, Year, `100mPlot`, Month_fct, n_samp_rnd) %>%
  mutate(index = row_number()) %>%
  left_join(A.estimates, by = "index") %>%
  filter(index %in% iter_idx) %>%
  ungroup() %>%
  group_by(Site, Year, Month_fct, n_samp_rnd, SiteTypeCode, index) %>%
  summarise(
    EffectiveAreaSampled = sum(A),
    PropAreaSampled = EffectiveAreaSampled / mean(TotalSiteArea),
    EffectiveNSampled = sum(value),
    EstTotalN = EffectiveNSampled / PropAreaSampled
  ) %>%
  ungroup() %>%
  group_by(Site, Year, Month_fct, n_samp_rnd, index) %>%
  summarise(EstTotalN = sum(EstTotalN)) %>%
  ungroup() %>%
  group_by(Site, Year, Month_fct) %>% # have summarized to month-level, but could incl n_samp_rnd here as well
  summarise(
    mean_N = mean(EstTotalN),
    sd_N = sd(EstTotalN),
    CV_N = sd_N / mean_N,
    l_95_N = quantile(EstTotalN, 0.025),
    l_90_N = quantile(EstTotalN, 0.05),
    l_50_N = quantile(EstTotalN, 0.25),
    med_N = quantile(EstTotalN, 0.5),
    u_50_N = quantile(EstTotalN, 0.75),
    u_90_N = quantile(EstTotalN, 0.95),
    u_95_N = quantile(EstTotalN, 0.975)
  )

write.csv(AbundanceEstimates_3, here("Results", "sqrt", "reanalysis_predictions_short_sqrt_upscaled.csv"))

AbundanceEstimates_3 <- AbundanceEstimates_3 %>%
  mutate(my = as.Date(paste(Year, as.numeric(Month_fct) + 3, "01", sep = "-")))

ggplot(AbundanceEstimates_3, aes(x = my, y = mean_N, colour = Site)) +
  geom_pointrange(aes(ymin = l_90_N, ymax = u_90_N)) +
  theme_bw() +
  ylab("Abundance estimate") +
  xlab("") +
  scale_color_discrete(name = "Population") +
  facet_wrap(~Site, ncol = 1, scales = "free") +
  theme(legend.position = "none")
ggsave(here("Results", "sqrt", "Scaled abundance estimates SHORT.jpg"),
  width = 1400, height = 1600, units = "px", dpi = 144
)

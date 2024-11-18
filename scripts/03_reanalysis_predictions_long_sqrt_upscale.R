#--------------------------------------------------------------------------- ---
# Code to upscale abundance, long timescale
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

load(file.path(repo_dir, "Data", "DH_data_all_4.rda"))

SiYr_all <- dh_dat_all_PPM_4 %>%
  select(Site, Year, `100mPlot`) %>%
  distinct() %>%
  mutate(across(-Site, factor))
SiYr <- nrow(SiYr_all)

dh <- suppressWarnings(
  pblapply(1:SiYr, function(ii) {
    temp <- dh_dat_all_PPM_4 %>%
      filter(
        Year == SiYr_all$Year[ii],
        Site == SiYr_all$Site[ii],
        `100mPlot` == SiYr_all$`100mPlot`[ii]
      ) %>%
      select(SubPlot, PPM_1:PPM_8) %>%
      filter(if_any(-1, ~ !is.na(.))) %>% # added here since some subplots in a grid not sampled in all years
      as.data.frame()
    temp <- Filter(function(x) !all(is.na(x)), temp)
    return(temp)
  })
)
names(dh) <- paste(SiYr_all$Year, SiYr_all$Site, SiYr_all$`100mPlot`, sep = "_")

k <- pmax(0, sapply(dh, ncol) - 1)

s <- sapply(dh, nrow)

# n_det[ii] = total number of detections
# s_dot[ii] = number of subplots with >0 detection
s_dot <- sapply(dh, function(cc) sum(rowSums(as.data.frame(cc[, -1])) > 0)) ## first column is SubPlot name
n_det <- sapply(dh, function(cc) sum(cc[, -1]))

sel <- s > 0 & k > 0 # AEB, changed here - filtering out rows where no survey done

site <- factor(SiYr_all$Site)
site_ref <- data.frame(level = levels(site), numeric = 1:length(levels(site)))
site <- as.numeric(site)

plot <- factor(SiYr_all$`100mPlot`)
plot_ref <- data.frame(level = levels(plot), numeric = 1:length(levels(plot)))
plot <- as.numeric(plot)

year <- factor(SiYr_all$Year)
year_ref <- data.frame(level = levels(year), numeric = 1:length(levels(year)))
year <- as.numeric(year)

SiYr_all <- SiYr_all %>%
  mutate(Year = as.numeric(as.character(Year)), `100mPlot` = as.character(`100mPlot`))
SiYr_wData <- SiYr_all[sel, ]
SiYr_wData <- SiYr_wData %>% left_join(X %>% select(Year, Site, `100mPlot`, SiteType, SiteTypeCode))

### Model results ----

load(here("reanalysis_long_samples_sqrt.Rdata"))

A.estimates <- do.call(rbind, post_subset(jagsfit, "A")) %>%
  as.data.frame() %>%
  mutate(index = row_number())

X %>%
  select(Site, `100mPlot`, SiteTypeCode, SiteType) %>%
  distinct() %>%
  group_by(Site, SiteTypeCode, SiteType) %>%
  summarise(n())

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
  group_by(Site, Year, `100mPlot`) %>%
  mutate(index = row_number()) %>%
  left_join(A.estimates, by = "index") %>%
  filter(index %in% iter_idx) %>%
  ungroup() %>%
  group_by(Site, Year, SiteTypeCode, index) %>%
  summarise(
    EffectiveAreaSampled = sum(A),
    PropAreaSampled = EffectiveAreaSampled / mean(TotalSiteArea),
    EffectiveNSampled = sum(value),
    EstTotalN = EffectiveNSampled / PropAreaSampled
  ) %>%
  ungroup() %>%
  group_by(Site, Year, index) %>%
  summarise(EstTotalN = sum(EstTotalN)) %>%
  ungroup() %>%
  group_by(Site, Year) %>%
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
AbundanceEstimates_3 %>% data.frame()

ggplot(AbundanceEstimates_3, aes(x = Year, y = mean_N, colour = Site)) +
  geom_pointrange(aes(ymin = l_90_N, ymax = u_90_N)) +
  theme_bw() +
  ylab("Abundance estimate") +
  scale_color_discrete(name = "Population") +
  scale_x_continuous(breaks = seq(2012, 2022, 2))
ggsave(here("Results", "sqrt", "Scaled abundance estimates.jpg"),
  width = 1200, height = 900, units = "px", dpi = 144
)

write.csv(AbundanceEstimates_3, here("Results", "sqrt", "reanalysis_predictions_sqrt_upscaled.csv"))

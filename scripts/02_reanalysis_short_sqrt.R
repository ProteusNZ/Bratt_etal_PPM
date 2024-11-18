#--------------------------------------------------------------------------- ---
# Code to run model, short timescale
# Client: Proteus
# Year: 2024
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

## LOAD DATA SETS ----

# load(here("Data", "DH_data_all.rda"))
# load(here("Data", "CR_tab.rda"))
# load(here("Data", "trap_rnd_dates.RData"))

## PREP OCCUPANCY DATA ----

# SiYr_all <- dh_dat_all_PPM %>%
#   select(Site, Year, `100mPlot`, Month_fct) %>%
#   distinct() %>% 
#   left_join(YSP %>% mutate(Site = toupper(Site), Month_fct = month-3), 
#             by = c("Site", "100mPlot", "Year", "Month_fct"),
#             ) %>% 
#   mutate(month = Month_fct + 3, 
#          n_samp_rnd = if_else(is.na(n_samp_rnd), 1, as.numeric(n_samp_rnd))) %>% 
#   uncount(n_samp_rnd) %>% 
#   group_by(Site, Year, `100mPlot`, Month_fct) %>% 
#   mutate(n_samp_rnd = row_number()) %>% 
#   ungroup() %>% 
#   mutate(across(-Site, factor))
# SiYr <- nrow(SiYr_all)
# 
# dh<-suppressWarnings(
#   pblapply(1:SiYr,function(ii){
#     temp<-dh_dat_all_PPM %>% filter(Year==SiYr_all$Year[ii],
#                                     Site==SiYr_all$Site[ii],
#                                     `100mPlot`==SiYr_all$`100mPlot`[ii],
#                                     Month_fct == SiYr_all$Month_fct[ii]
#                                     ) %>%
#       select(SubPlot, PPM_1:PPM_5 ) %>%
#       filter(if_any(-1, ~!is.na(.))) %>% # added here since some subplots in a grid not sampled in all years
#       as.data.frame()
#     temp<-Filter(function(x) !all(is.na(x)),temp) 
#     return(temp)
#   })
# )
# names(dh)<-paste(SiYr_all$Year,SiYr_all$Site,SiYr_all$`100mPlot`,SiYr_all$Month_fct,SiYr_all$n_samp_rnd,sep="_")
# 
# k<-pmax(0,sapply(dh,ncol)-1)
# 
# s<-sapply(dh,nrow)
# 
# # n_det[ii] = total number of detections
# # s_dot[ii] = number of subplots with >0 detection
# s_dot<-sapply(dh,function(cc) sum(rowSums(as.data.frame(cc[,-1]))>0)) ## first column is SubPlot name
# n_det<-sapply(dh,function(cc) sum(cc[,-1]))
# 
# sel<- s>0 & k > 0 # AEB, changed here - filtering out rows where no survey done
# 
# site<-factor(SiYr_all$Site)
# site_ref<-data.frame(level=levels(site),numeric=1:length(levels(site)))
# site<-as.numeric(site)
# 
# plot<-factor(SiYr_all$`100mPlot`)
# plot_ref<-data.frame(level=levels(plot),numeric=1:length(levels(plot)))
# plot<-as.numeric(plot)
# 
# year<-factor(SiYr_all$Year)
# year_ref<-data.frame(level=levels(year),numeric=1:length(levels(year)))
# year<-as.numeric(year)
# 
# month<-factor(SiYr_all$Month_fct)
# mont_ref<-data.frame(level=levels(month),numeric=1:length(levels(month)))
# month<-as.numeric(month)

## PREP CAPTURE DATA ----

# YSP2<-YSP %>% 
#   uncount(n_samp_rnd) %>% 
#   group_by(Site, Year, `100mPlot`, month) %>% 
#   mutate(n_samp_rnd = row_number()) %>% 
#   ungroup() %>% 
#   mutate(Site = toupper(Site), Month_fct = month-3)
# 
# # trap_rnd_dates<-trap_rnd_dates %>% 
# #   filter(month(mid_date)>=4,month(mid_date)<=7) %>% # AEB note - changed to midpoint here
# #   mutate(month = month(mid_date)) %>% # AEB note added month column here
# #   group_by(Site, `100mPlot`, Year, month) %>% 
# #   filter(row_number() == 1)
# 
# # CR_tab2<-CR_tab %>% 
# #   filter(Fate=="released") %>% 
# #   #semi_join(trap_rnd_dates,by=c("Site","100mPlot","Year","samp_rnd")) %>% 
# #   # this achieves the same join, but including the month variable
# #   semi_join(trap_rnd_dates %>% select(Site, `100mPlot`, Year, samp_rnd, month), 
# #             by=c("Site","100mPlot","Year","samp_rnd", "month")) %>% 
# #   filter(!is.na(month)) %>% 
# #   mutate(caps = as.numeric(rowSums(data.frame(A1,A2,A3,B1,B2,B3),na.rm=TRUE)>0)) %>%
# #   group_by(Site, `100mPlot`, Year, Unique_Indiv_ID, month) %>% 
# #   mutate(caps2 = as.numeric(sum(caps)>0)) %>% 
# #   distinct(Site, `100mPlot`, Year, Unique_Indiv_ID, month, caps2) %>% 
# #   arrange(month, Site, `100mPlot`, Year, Unique_Indiv_ID) %>% 
# #   pivot_wider(id_cols = c("Site","100mPlot","Year", "Unique_Indiv_ID"),
# #              names_from = month,values_from = caps2, values_fill = 0) %>%
# #   pivot_longer(cols=as.character(4:7),names_to = "month",values_to = "caps") %>%
# #   mutate(month=as.numeric(month)-3) %>%
# #   left_join(trap_rnd_dates %>% select(Site,`100mPlot`,Year,month,min_date) %>% mutate(month = month-3)) %>%
# #   mutate(caps=replace(caps,is.na(min_date),NA))
# 
# # CR_tab2 <- CR_tab2 %>%
# #   pivot_wider(id_cols = c("Site","100mPlot","Year","Unique_Indiv_ID"),
# #               names_from = month,values_from = caps) %>% 
# #   mutate(Site = toupper(Site))
# 
# CR_tab2 <- CR_tab %>% 
#   mutate(Site = toupper(Site), 
#          Month_fct = month-3)
# 
# sel_sess<-paste0(rep(LETTERS[1:2],each=3),rep(1:3,2))
# 
# # AEB note - now loops over SiYr, so occupancy and capture indices match
# ch<- pblapply(1:nrow(SiYr_all),function(ii){ 
#   # lapply(1:nrow(YSP2),function(ii){ 
#   temp <- CR_tab2 %>% filter(Year==SiYr_all$Year[ii],
#                              Site==SiYr_all$Site[ii],
#                              `100mPlot`==SiYr_all$`100mPlot`[ii],
#                              Month_fct == SiYr_all$Month_fct[ii], 
#                              samp_rnd == SiYr_all$n_samp_rnd[ii]) %>%
#     select(Unique_Indiv_ID, all_of(sel_sess)) %>% 
#     distinct() %>% 
#     as.data.frame()
#   # temp<-Filter(function(x) !all(is.na(x)),temp)
#   return(temp)
# })
# names(ch)<-paste(SiYr_all$Year,SiYr_all$Site,SiYr_all$`100mPlot`,SiYr_all$Month_fct,sep="_")
# 
# sum(sapply(ch,nrow))
# 
# ## summary stats
# # t <-YSP %>% filter(!is.na(samp_rnd)) %>% select(!month) %>%
# #   semi_join(trap_rnd_dates) %>%
# #   pivot_wider(names_from = samp_rnd,values_from = t, values_fill = NA)
# 
# t<-SiYr_all %>% 
#   mutate(Year = as.numeric(as.character(Year)), 
#          `100mPlot` = as.character(`100mPlot`),
#          month = as.numeric(as.character(month)), 
#          n_samp_rnd = as.numeric(as.character(n_samp_rnd))) %>% 
#   ungroup() %>% 
#   rowwise() %>% 
#   mutate(inYSP = if_else(paste(Site, `100mPlot`, Year, month, n_samp_rnd, sep = "-") %in% do.call(paste, c(YSP2[1:5], sep="-")), 
#                          1, 0)) %>% 
#   mutate(t = if_else(inYSP == 1, 6, 0)) # works since no NAs at this scale
# 
# t<-t$t
# 
# Mt1<-sapply(ch,nrow)
# n_cap<-(sapply(ch,function(cc) sum(cc[,-1]))) ## first column is id
# cbind(Mt1,n_cap,t)
# 
# n_cap<-replace(n_cap,is.na(n_cap),0)
# 
# t_max <- 6 
# 
# save.image(here("Data","short_pre_jags_image.RData"))

## LOAD PREPROCESSED DATA ----

load(here("Data", "short_pre_jags_image.RData"))
load(here("Data", "corrected_X.rda"))

## SET UP MODEL ----

# SiYr = number of site/year data sets
# t[ii] = number of capture occasions
# n_cap[ii] = number of captures
# Mt1[ii] = number of unique indviduals
# hug_ones[ii] = vector of 1's for the ones-trick
# n_det[ii] = total number of detections
# s_dot[ii] = number of subplots with >0 detection
# k[ii] = number of detection occasions 
# s[ii] = number of subplots
# occ_ones[ii] = vector of 1's for the ones-trick
mcmc.data<-list(SiYr = sum(sel),
                

                t = t[sel],
                #t_max = t_max,
                n_cap = n_cap[sel],
                Mt1 = Mt1[sel],
                hug_ones = rep(1,SiYr)[sel],
                
                n_det = n_det[sel],
                s_dot = s_dot[sel],
                k = k[sel],
                s = s[sel],
                occ_ones = rep(1,SiYr)[sel],
                
                plot=plot[sel],n_plot=max(plot),
                month=month[sel],n_month=max(month),
                year=year[sel],n_year=max(year),

                site_type = X$SiteTypeCode[sel],
                d=11.52,se_d=1.04,L=37.5,pi=pi  # 1 month
                
)

mcmc.pars<-c("beta","sd_eps", 
             "hug_p_sd","occ_p_sd","psi_sd",
             "mu_hug_p", "mu_occ_p", "mu_psi",
             "plot_occ_p_sd","plot_occ_p",
             "month_occ_p_sd","month_occ_p",
             "year_occ_p_sd","year_occ_p",
             "plot_hug_p_sd","plot_hug_p",
             "month_hug_p_sd","month_hug_p",
             "year_hug_p_sd","year_hug_p", 
             "plot_psi_sd","plot_psi",
             "month_psi_sd","month_psi",
             "year_psi_sd","year_psi",
             "beta_type_occ_p", "beta_type_psi",
             "den","A",
             "hug_p","occ_p","psi",
             "prop_occ",
             "N"
)

jags.mod <- jags.model(file.path(repo_dir,"Code", "Reanalysis" ,"reanalysis_short_jags_sqrt.R"), mcmc.data, n.chains = 3)

## RUN MODEL ----

t.start <- Sys.time()
update(jags.mod,10000)
t.end <- Sys.time()
(runTime <- t.end - t.start)

t.start <- Sys.time()
jagsfit <- coda.samples(jags.mod, variable.names = mcmc.pars, n.iter=300000, thin=10)
t.end <- Sys.time()
(runTime <- t.end - t.start)
beep(sound = 8)

## PROCESS AND SAVE RESULTS ----

save(jagsfit, file = file.path(SP_dir,"Results","short","reanalysis_short_samples_sqrt.Rdata"))

# get summary stats
summ <- t(post_summ(jagsfit, get_params(jagsfit, type = "base_index"), 
                    neff = TRUE, Rhat = TRUE, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))) %>% 
  as.data.frame()
summ <- summ %>% 
  rownames_to_column(var = "name")
save(summ, file = file.path(SP_dir,"Results","short","reanalysis_short_summary_sqrt.Rdata"))

model.pars <- post_subset(jagsfit, c("beta\\[", "sd_eps"))
hug_p.pars <- post_subset(jagsfit, c("mu_hug_p", "hug_p_sd", 
                                     "^plot_hug_p", "plot_hug_p_sd", 
                                     "^year_hug_p", "year_hug_p_sd", 
                                     "^hug_p"))
occ_p.pars <- post_subset(jagsfit, c("mu_occ_p", "occ_p_sd", 
                                     "^plot_occ_p", "plot_occ_p_sd", 
                                     "^year_occ_p", "year_occ_p_sd", 
                                     "^month_occ_p", "month_occ_p_sd", 
                                     "^occ_p"))
psi.pars <- post_subset(jagsfit, c("mu_psi", "psi_sd", 
                                   "^plot_psi", "plot_psi_sd", 
                                   "^year_psi", "year_psi_sd", 
                                   "^month_psi", "month_psi_sd", 
                                   "^psi")) 
cov.pars <- post_subset(jagsfit, c("beta_"))

save(model.pars, file=file.path(SP_dir,"Results","short","reanalysis_short_model_pars_sqrt.RData"))
save(hug_p.pars, file = file.path(SP_dir,"Results","short","reanalysis_short_hug_p_pars_sqrt.RData"))
save(occ_p.pars, file = file.path(SP_dir,"Results","short","reanalysis_short_occ_p_pars_sqrt.RData"))
save(psi.pars, file = file.path(SP_dir,"Results","short","reanalysis_short_psi_pars_sqrt.RData"))
save(cov.pars, file = file.path(SP_dir,"Results","short","reanalysis_short_cov_pars_sqrt.RData"))

# the 501 combos are SiYr_all[sel, ]
# the 70 we care about right now are YSP2
# so do some indexing magic
est <- summ %>% 
  filter(!str_detect(name, "sd") &
           (str_detect(name, "den") |
              str_detect(name, "N") |
              str_detect(name, "^hug_p") |
              str_detect(name, "^occ_p") |
              str_detect(name, "^psi") |
              str_detect(name, "prop_occ"))
  ) 
colsToUse <- intersect(colnames(SiYr_all), colnames(YSP2))
matches <- match(do.call("paste", YSP2[, colsToUse]), do.call("paste", SiYr_all[sel, colsToUse]))
matches <- bind_cols(YSP = 1:length(matches), index = matches) %>% 
  filter(!is.na(index)) %>%
  as.data.frame()
n_cap2 <- n_cap[sel]
t2 <- t[sel]
est <- est %>% 
  mutate(index = str_first_number(name)) %>% 
  filter(index %in% matches$index) %>% 
  left_join(matches, by = "index") %>% 
  left_join(YSP2 %>% 
              rownames_to_column(var = "YSP") %>% 
              mutate(YSP = as.numeric(YSP)), by = "YSP") %>% 
  rename(l_95 = `2.5%`, l_90 = `5%`, l_50 = `25%`, med = `50%`, u_50 = `75%`, u_90 = `95%`, u_95 = `97.5%`) %>% 
  mutate(name = gsub("\\s*\\[[^\\)]+\\]","",name)) %>% 
  select(Site, `100mPlot`, Year, name,month,n_samp_rnd,Month_fct, 2:14) %>% 
  pivot_wider(id_cols = c(Site, `100mPlot`, Year,month,n_samp_rnd,Month_fct), 
              names_from = name, 
              values_from = 8:16, 
              names_glue = "{name}_{.value}", 
              names_vary = "slowest") %>% 
  bind_cols(s_dot = s_dot[sel][matches$index], 
            n_det = n_det[sel][matches$index], 
            Mt1 = Mt1[sel][matches$index], 
            n_cap_total = n_cap2[matches$index], 
            t_total = t2[matches$index]) 

write.csv(est,file.path(SP_dir,"Results","short","short time density and occupancy estimates_reanalysis_sqrt.csv"))



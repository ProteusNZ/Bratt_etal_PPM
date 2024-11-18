# JAGS model for basic Huggins CR model and Occupancy model ####

model{
  
  beta[1] ~ dnorm(-3,1) # intercept
  beta[2] ~ dnorm(0,1/4) # prop_occ effect
  beta[3] ~ dnorm(0,1/4) # occ_p effect
  
  
  sd_eps ~ dnorm(0, 1/0.5625)T(0, 1.5)
  tau_eps <- 1/(sd_eps^2)
  
  hug_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  hug_p_tau <- 1/(pow(hug_p_sd,2))
  
  occ_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  occ_p_tau <- 1/(pow(occ_p_sd,2))
  
  psi_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  psi_tau <- 1/(pow(psi_sd,2))
  
  ## means for abund, detection, occ
  mu_hug_p ~ dnorm(-1, 1)
  mu_occ_p ~ dnorm(0, 1)
  mu_psi ~ dnorm(0, 1)
  
  ## plot effects
  plot_occ_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  plot_occ_p_tau <- 1/(pow(plot_occ_p_sd,2))
  plot_hug_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  plot_hug_p_tau <- 1/(pow(plot_hug_p_sd,2))
  plot_psi_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  plot_psi_tau <- 1/(pow(plot_psi_sd,2))
  for(ii in 1:n_plot){
    plot_occ_p[ii] ~ dnorm(0,plot_occ_p_tau)
    plot_hug_p[ii] ~ dnorm(0,plot_hug_p_tau)
    plot_psi[ii] ~ dnorm(0,plot_psi_tau)
  }
  
  
  ## month effects
  month_occ_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  month_occ_p_tau <- 1/(pow(month_occ_p_sd,2))
  month_hug_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  month_hug_p_tau <- 1/(pow(month_hug_p_sd,2))
  month_psi_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  month_psi_tau <- 1/(pow(month_psi_sd,2))
  for(ii in 1:n_month){
    month_occ_p[ii] ~ dnorm(0,month_occ_p_tau)
    month_hug_p[ii] ~ dnorm(0,month_hug_p_tau)
    month_psi[ii] ~ dnorm(0,month_psi_tau)
  }
  
  
  ## year effects
  year_occ_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  year_occ_p_tau <- 1/(pow(year_occ_p_sd,2))
  year_hug_p_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  year_hug_p_tau <- 1/(pow(year_hug_p_sd,2))
  year_psi_sd ~ dnorm(0, 1/0.5625)T(0, 1.5)
  year_psi_tau <- 1/(pow(year_psi_sd,2))
  for(ii in 1:n_year){
    year_occ_p[ii] ~ dnorm(0,year_occ_p_tau)
    year_hug_p[ii] ~ dnorm(0,year_hug_p_tau)
    year_psi[ii] ~ dnorm(0,year_psi_tau)
  }
  
  beta_type_occ_p ~ dnorm(0, 1/0.25)
  beta_type_psi ~ dnorm(0, 1/0.25)
  
  ## end priors
  
  ## area calculation
  inv_var_se_d <- 1/(se_d*se_d)
  d_star ~ dnorm(d,inv_var_se_d)T(0,)
  W<-d_star/2
  A<-(L*L+4*L*W+pi*W*W)/10000
  
  ## likelihood
  for(ii in 1:SiYr){
    eps[ii] ~ dnorm(0,tau_eps)
    
    log(den[ii]) <- beta[1] + beta[2]*sqrt(prop_occ[ii]) + beta[3]*sqrt(occ_p[ii]) + eps[ii]
    
    hug_p_mu[ii] <- mu_hug_p + plot_hug_p[plot[ii]]+ month_hug_p[month[ii]] + year_hug_p[year[ii]]
    
    ### note change for site type
    occ_p_mu[ii] <- mu_occ_p + beta_type_occ_p*(site_type[ii]-1) + plot_occ_p[plot[ii]] + month_occ_p[month[ii]] + year_occ_p[year[ii]]
    psi_mu[ii] <- mu_psi + beta_type_psi*(site_type[ii]-1) + plot_psi[plot[ii]] + month_psi[month[ii]] + year_psi[year[ii]]
    
    l_hp[ii] ~ dnorm(hug_p_mu[ii],hug_p_tau)
    logit(hug_p[ii]) <- l_hp[ii]
    
    l_op[ii] ~ dnorm(occ_p_mu[ii],occ_p_tau)
    logit(occ_p[ii]) <- l_op[ii]
    
    l_psi[ii] ~ dnorm(psi_mu[ii],psi_tau)
    logit(psi[ii]) <- l_psi[ii]
    
    ## occupancy model
    occ_like[ii] <- psi[ii]^s_dot[ii]*occ_p[ii]^(n_det[ii])*(1-occ_p[ii])^(s_dot[ii]*k[ii]-n_det[ii]) * (1 - psi[ii]*(1-(1-occ_p[ii])^k[ii]))^(s[ii]-s_dot[ii])
    occ_ones[ii]~dbern(occ_like[ii])
    
    psi_con[ii] <- psi[ii]*((1-occ_p[ii])^k[ii])/(1-psi[ii]*(1-(1-occ_p[ii])^k[ii]))
    s_0[ii] ~ dbin(psi_con[ii],s[ii]-s_dot[ii])
    prop_occ[ii] <- (s_dot[ii]+s_0[ii])/s[ii]
    
    ## Huggins model
    ### full likelihood
    p_star[ii]<-1-(1-hug_p[ii])^t[ii]
    lambda[ii] <- den[ii]*A
    
    f0[ii] ~ dpois(lambda[ii]*(1-p_star[ii]))
    Mt1[ii] ~ dpois(lambda[ii]*p_star[ii])
    N[ii] <- Mt1[ii] + f0[ii]
    
    full_like[ii]<-hug_p[ii]^(n_cap[ii])*(1-hug_p[ii])^(Mt1[ii]*t[ii] - n_cap[ii])/(p_star[ii]^Mt1[ii])
    hug_ones[ii]~dbern(full_like[ii])
    
  }
  
}
################################################################################################################
################################################################################################################
rm(list=ls())
library(base)
library(boot)
library(bayesplot)
library(posterior)
library(rjags)
library(coda)
library(R2jags)
library(jagstools)
library(mcmcplots)
library(tidylog)
library(tidyverse)
################################################################################################################
################################################################################################################
################################################################################################################

options(mc.cores = parallel::detectCores())

#set your working directory
data <- data.frame(read_csv("example_data.csv"))
colnames(data)
no_of_years <- 9
no_counts_perFY <- length(which(data$FinancialYear=='2014/15'))

data_list <- list(
N_size=dim(data)[1],
age_g1=c(ifelse(data$AgeGroup=='15-34',1,0)),
age_g2=c(ifelse(data$AgeGroup=='35-49',1,0)),
age_g3=c(ifelse(data$AgeGroup=='50-64',1,0)),
sex_g=c(ifelse(data$Sex=='Males',1,0)),
year2014=c(ifelse(data$FinancialYear=='2014/15',1,0)),
year2015=c(ifelse(data$FinancialYear=='2015/16',1,0)),
year2016=c(ifelse(data$FinancialYear=='2016/17',1,0)),
year2017=c(ifelse(data$FinancialYear=='2017/18',1,0)),
year2018=c(ifelse(data$FinancialYear=='2018/19',1,0)),
year2019=c(ifelse(data$FinancialYear=='2019/20',1,0)),
year2020=c(ifelse(data$FinancialYear=='2020/21',1,0)),
year2021=c(ifelse(data$FinancialYear=='2021/22',1,0)),
year2022=c(ifelse(data$FinancialYear=='2022/23',1,0)),
year=c(rep(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9),3)),
P=c(data$Population),
n_known=c(data$NumberBaselineCohort),
n_known_on=c(data$NumberBaselineCohortOnOat),
pyr_on=c(data$PersonYearsAtRiskCohortOnOat),
pyr_off=c(data$PersonYearsAtRiskCohortOffOat),
events_linked_on=cbind(c(data$DeathsCohortOnOat),c(data$HospitalisationsCohortOnOat)),
events_linked_off=cbind(c(data$DeathsCohortOffOat),c(data$HospitalisationsCohortOffOat)),
events_extra=cbind(c(data$DeathsUnobserved),c(data$HospitalisationsUnobserved)),
other_cause_exit=c(data$OtherCauseMortality),
pyr_died=c(data$PersonYearsAtRiskDiedDrd)
)
########################################################################
########################################################################
### Initial Values

prev_lower_bound <- data_list$n_known/data_list$P
n_main_parameters <- 12
N_size <- data_list$N_size
baseline_location <- which(data_list$year2014==1 & data_list$sex==0 & data_list$age_g2==0 & data_list$age_g3==0)






initial_values <- list(
  list(
    beta0=rep(-4.20,2),
    beta=cbind(c(0,rep(0.5,n_main_parameters-1)),c(0,rep(0.5,n_main_parameters-1))),
    beta_h=rep(0,3),
    EvRate_RE_treat_fy = matrix(0,nrow=no_of_years,ncol=2),
    sd_ER_re_treat_fy = rep(1.5,2),
    #####################################
    beta_others=c(1.0,rep(0,n_main_parameters-1)),
    #####################################
#    prev=(prev_lower_bound+0.0001)*1.2,
    gamma=c(logit((prev_lower_bound+0.0001)*1.2)[baseline_location],rep(0,n_main_parameters-1)),
    gamma_agesex=rep(0,2),
    sd_PREV_RE_age2_fy=1.5,
    sd_PREV_RE_age3_fy=1.5,
    PREV_RE_age2_fy=rep(0,no_of_years),
    PREV_RE_age3_fy=rep(0,no_of_years),
    sd_PREV_RE_sex_fy=1.5,
    PREV_RE_sex_fy=rep(0,no_of_years),
    prev_known=prev_lower_bound+0.0001
  ) 
  ,  
  list(
    beta0=rep(-4.05,2),
    beta=cbind(c(0,rep(0,n_main_parameters-1)),c(0,rep(0,n_main_parameters-1))),
    beta_h=rep(0,3),
    EvRate_RE_treat_fy = matrix(0,nrow=no_of_years,ncol=2),
    sd_ER_re_treat_fy = rep(2.5,2),
    #####################################
    beta_others=c(0.8,rep(1,n_main_parameters-1)),
    #####################################
#    prev=(prev_lower_bound+0.0001)*1.5,
    gamma=c(logit((prev_lower_bound+0.0001)*1.5)[baseline_location],rep(1,n_main_parameters-1)),
    gamma_agesex=rep(0,2),
    sd_PREV_RE_age2_fy=2.5,
    sd_PREV_RE_age3_fy=2.5,
    PREV_RE_age2_fy=rep(0,no_of_years),
    PREV_RE_age3_fy=rep(0,no_of_years),
    sd_PREV_RE_sex_fy=2.5,
    PREV_RE_sex_fy=rep(0,no_of_years),
    prev_known=prev_lower_bound+0.0001
  )
  ,  
  list(
    beta0=rep(-4.35,2),
    beta=cbind(c(0,rep(-0.5,n_main_parameters-1)),c(0,rep(-0.5,n_main_parameters-1))),
    beta_h=rep(0,3),
    EvRate_RE_treat_fy = matrix(0,nrow=no_of_years,ncol=2),
    sd_ER_re_treat_fy = rep(1.9,2),
    #####################################
    beta_others=c(0.6,rep(0.1,n_main_parameters-1)),
    #####################################
#    prev=(prev_lower_bound+0.0001)*2,
    gamma=c(logit((prev_lower_bound+0.0001)*2)[baseline_location],rep(0.1,n_main_parameters-1)),
    gamma_agesex=rep(0,2),
    sd_PREV_RE_age2_fy=1.9,
    sd_PREV_RE_age3_fy=1.9,
    PREV_RE_age2_fy=rep(0,no_of_years),
    PREV_RE_age3_fy=rep(0,no_of_years),
    sd_PREV_RE_sex_fy=1.9,
    PREV_RE_sex_fy=rep(0,no_of_years),
    prev_known=prev_lower_bound+0.0001
  )
)


########################################################################
########################################################################


jags_model <- function(){
#####################
#  MORTALITY MODEL  #
#####################
for(i in 1:N_size){
  for(event in 1:2){
  ## Likelihood:
  events_linked_off[i,event] ~ dpois(mean_events_linked_off[i,event])
  events_linked_on[i,event] ~ dpois(mean_events_linked_on[i,event])  
  
  mean_events_linked_off[i,event] <- lambda_off[i,event] * pyr_off[i]
  mean_events_linked_on[i,event] <- lambda_on[i,event] * pyr_on[i]
  
  log(lambda_off[i,event]) <- beta0[event] + 
              beta[2,event]*age_g2[i] +     
              beta[3,event]*age_g3[i] +
              beta[4,event]*sex_g[i] +
              beta[5,event]*year2015[i] +
              beta[6,event]*year2016[i] +
              beta[7,event]*year2017[i] +
              beta[8,event]*year2018[i] +
              beta[9,event]*year2019[i] +
              beta[10,event]*year2020[i] +
              beta[11,event]*year2021[i] +
              beta[12,event]*year2022[i] +
              beta_h[1]*age_g2[i]*sex_g[i]*(event-1) +
              beta_h[2]*age_g3[i]*sex_g[i]*(event-1)
  log(lambda_on[i,event]) <- log(lambda_off[i,event]) +
                  beta[1,event] +
                  beta_h[3]*sex_g[i]*(event-1) +
                  EvRate_RE_treat_fy[year[i],event]
   }
}

## Priors for mortality model:
## intercept 
for(event in 1:2){
  beta0[event] ~ dnorm(-4.1, 1/100)
}
## Main effects age, gender, treatment
for(j in 1:12){
  for(event in 1:2){
    beta[j,event] ~ dnorm(0, 1/100)
  }
}
beta_h[1] ~ dnorm(0,1/100)
beta_h[2] ~ dnorm(0,1/100)
beta_h[3] ~ dnorm(0,1/100)
for(y in 1:9){
    for(event in 1:2){
      EvRate_RE_treat_fy[y,event] ~ dnorm(0,pow(sd_ER_re_treat_fy[event],-2))
    }
  }
  sd_ER_re_treat_fy[1] ~ dunif(0,5)
  sd_ER_re_treat_fy[2] ~ dunif(0,5)
#######################################################
##  MORTALITY RATE MODEL # Other Deaths, ACM - HMBBs  #
#######################################################
for(i in 1:N_size){
  other_cause_exit[i] ~ dpois(mean_other_death[i])
  mean_other_death[i] <- lambda_others[i]*pyr_off[i]
  log(lambda_others[i]) <- beta_others[1] + 
               beta_others[2]*age_g2[i] +     
               beta_others[3]*age_g3[i] +
               beta_others[4]*sex_g[i] +
               beta_others[5]*year2015[i] +
               beta_others[6]*year2016[i] +
               beta_others[7]*year2017[i] +
               beta_others[8]*year2018[i] +
               beta_others[9]*year2019[i] +
               beta_others[10]*year2020[i] +
               beta_others[11]*year2021[i] +
               beta_others[12]*year2022[i]
}

## Priors for mortality model:
## intercept
beta_others[1] ~ dnorm(1, 1/100)
## Main effects age, gender, treatment
for(j in 2:12){
  beta_others[j] ~ dnorm(0, 1/100)
}

######################
# PREVALENCE MODEL   #
######################
  for(i in 1:N_size){
    events_extra[i,1] ~ dpois(mean_events_extra[i,1]) 
    mean_events_extra[i,1] <- lambda_off[i,1] * pyr_extra[i] 
    events_extra[i,2] ~ dpois(mean_events_extra[i,2]) 
    mean_events_extra[i,2] <- lambda_off[i,2] * pyr_extra[i] 
    correction_term[i] <- (1/lambda_others[i])*(1 - exp(-lambda_others[i]))
    pyr_extra[i] <- pyr_died[i] + (n_extra[i] - pyr_died[i])*correction_term[i]
    n_extra[i] <- N[i] - n_known[i]
    N[i] <- P[i] * prev[i]
    prev[i] <- ilogit(gamma[1] +
                      gamma[2]*age_g2[i] +
                      gamma[3]*age_g3[i] +
                      gamma[4]*sex_g[i] + 
                      gamma[5]*year2015[i] +
                      gamma[6]*year2016[i] +
                      gamma[7]*year2017[i] +
                      gamma[8]*year2018[i] +
                      gamma[9]*year2019[i] +
                      gamma[10]*year2020[i] +
                      gamma[11]*year2021[i] +
                      gamma[12]*year2022[i] +
                      PREV_RE_age2_fy[year[i]]*age_g2[i] +
                      PREV_RE_age3_fy[year[i]]*age_g3[i]) + 
                      prev_known[i] 
    n_known[i] ~ dbin(prev_known[i],P[i])
    prev_known[i] ~ dunif(0,1)
  }
    gamma[1] ~ dnorm(-4.2,1/10)
  for(i in 2:12){
    gamma[i] ~ dnorm(0,1/10)
  }
  for(y in 1:9){
    PREV_RE_age2_fy[y] ~ dnorm(0,pow(sd_PREV_RE_age2_fy,-2))
    PREV_RE_age3_fy[y] ~ dnorm(0,pow(sd_PREV_RE_age3_fy,-2))
  }
  sd_PREV_RE_age2_fy ~ dunif(0,3)
  sd_PREV_RE_age3_fy ~ dunif(0,3)
}




gc()
set.seed(123)
jagsfitted <- R2jags::jags(data = data_list, inits = initial_values,
   parameters.to.save =c("beta0","beta","beta_h","mean_events_linked_on","mean_events_linked_on",
                         "beta_others","mean_other_death","gamma","mean_events_extra","N","prev",
                         "sd_PREV_RE_age2_fy","sd_PREV_RE_age3_fy","sd_ER_re_treat_fy"), 
   n.chains = 3, n.iter=20000, n.burnin = 10000, n.thin=1,
   DIC=FALSE, refresh = 1000,
   model.file = jags_model)








jagsfit_table<-jagsfitted$BUGSoutput$summary
range(jagsfit_table[,"n.eff"])
range(jagsfit_table[,"Rhat"])


########################################################################
########################################################################
## draw prevalence estimates per Financial Year
jagsmcmc <- as.mcmc(jagsfitted)
posterior_sample <- data.frame(rbind(jagsmcmc[[1]],jagsmcmc[[2]],jagsmcmc[[3]]))

P <- data_list$P
N_size <- data_list$N_size
vector2014 <- c(which(data_list$year==1))
vector2015 <- c(which(data_list$year==2))
vector2016 <- c(which(data_list$year==3))
vector2017 <- c(which(data_list$year==4))
vector2018 <- c(which(data_list$year==5))
vector2019 <- c(which(data_list$year==6))
vector2020 <- c(which(data_list$year==7))
vector2021 <- c(which(data_list$year==8))
vector2022 <- c(which(data_list$year==9))


draws <- data.frame(posterior_sample[,paste("N.",1:N_size,".",sep="")])

draws <- draws %>%
    mutate(prev_2014 = apply(draws[,vector2014],1,sum)/sum(data_list$P[vector2014])*100) %>%
    mutate(prev_2015 = apply(draws[,vector2015],1,sum)/sum(data_list$P[vector2015])*100) %>%
    mutate(prev_2016 = apply(draws[,vector2016],1,sum)/sum(data_list$P[vector2016])*100) %>%
    mutate(prev_2017 = apply(draws[,vector2017],1,sum)/sum(data_list$P[vector2017])*100) %>%
    mutate(prev_2018 = apply(draws[,vector2018],1,sum)/sum(data_list$P[vector2018])*100) %>%
    mutate(prev_2019 = apply(draws[,vector2019],1,sum)/sum(data_list$P[vector2019])*100) %>%
    mutate(prev_2020 = apply(draws[,vector2020],1,sum)/sum(data_list$P[vector2020])*100) %>%
    mutate(prev_2021 = apply(draws[,vector2021],1,sum)/sum(data_list$P[vector2021])*100) %>%
    mutate(prev_2022 = apply(draws[,vector2022],1,sum)/sum(data_list$P[vector2022])*100) %>%
    tidylog::select(prev_2014,prev_2015,prev_2016,prev_2017,prev_2018,prev_2019,
                    prev_2020,prev_2021,prev_2022)

## summarise results
summarise_draws(draws, mean, ~quantile(.x, probs = c(0.025, 0.975)))



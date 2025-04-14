################################################################################################################
################################################################################################################
rm(list=ls())

library(cmdstanr)
library(rstan)
library(bayesplot)
library(loo)
library(posterior)
library(readr)
library(base)
library(boot)
library(tidyverse)
################################################################################################################
################################################################################################################
################################################################################################################


options(mc.cores = parallel::detectCores())

#set your working directory
stan_data <- data.frame(read_csv("example_data.csv"))
colnames(stan_data)
no_of_years <- 9
no_counts_perFY <- length(which(stan_data$FinancialYear=='2014/15'))

data_list <- list(
N_size=dim(stan_data)[1],
age_g1=c(ifelse(stan_data$AgeGroup=='15-34',1,0)),
age_g2=c(ifelse(stan_data$AgeGroup=='35-49',1,0)),
age_g3=c(ifelse(stan_data$AgeGroup=='50-64',1,0)),
sex_g=c(ifelse(stan_data$Sex=='Males',1,0)),
year2014=c(ifelse(stan_data$FinancialYear=='2014/15',1,0)),
year2015=c(ifelse(stan_data$FinancialYear=='2015/16',1,0)),
year2016=c(ifelse(stan_data$FinancialYear=='2016/17',1,0)),
year2017=c(ifelse(stan_data$FinancialYear=='2017/18',1,0)),
year2018=c(ifelse(stan_data$FinancialYear=='2018/19',1,0)),
year2019=c(ifelse(stan_data$FinancialYear=='2019/20',1,0)),
year2020=c(ifelse(stan_data$FinancialYear=='2020/21',1,0)),
year2021=c(ifelse(stan_data$FinancialYear=='2021/22',1,0)),
year2022=c(ifelse(stan_data$FinancialYear=='2022/23',1,0)),
year=c(rep(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9),3)),
P=c(stan_data$Population),
n_known=c(stan_data$NumberBaselineCohort),
n_known_on=c(stan_data$NumberBaselineCohortOnOat),
pyr_on=c(stan_data$PersonYearsAtRiskCohortOnOat),
pyr_off=c(stan_data$PersonYearsAtRiskCohortOffOat),
events_linked_on=cbind(c(stan_data$DeathsCohortOnOat),c(stan_data$HospitalisationsCohortOnOat)),
events_linked_off=cbind(c(stan_data$DeathsCohortOffOat),c(stan_data$HospitalisationsCohortOffOat)),
events_extra=cbind(c(stan_data$DeathsUnobserved),c(stan_data$HospitalisationsUnobserved)),
other_cause_exit=c(stan_data$OtherCauseMortality),
pyr_died=c(stan_data$PersonYearsAtRiskDiedDrd)
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
    prev=(prev_lower_bound+0.0001)*1.2,
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
    prev=(prev_lower_bound+0.0001)*1.5,
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
    prev=(prev_lower_bound+0.0001)*2,
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


stanmod <- cmdstan_model("model/mpep_model.stan", compile = TRUE)


stanfitted <- stanmod$sample(data=data_list,
                             seed=123,
                             iter_warmup=2000,
                             iter_sampling=2000,
                             chains = 3,
                             parallel_chains=3,
                             save_warmup = TRUE,
                             thin=1,
                             max_treedepth=15,
                             init=initial_values)
stanfit_table <- data.frame(stanfitted$summary())
range(na.omit(stanfit_table$rhat))
range(na.omit(stanfit_table$ess_tail))
range(na.omit(stanfit_table$ess_bulk))


########################################################################
########################################################################
## draw prevalence estimates per Financial Year
posterior_sample <- data.frame(stanfitted$draws(format = "matrix",inc_warmup = FALSE))

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



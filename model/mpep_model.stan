//opi_oid dependence model - with treatment effect

data{
  int<lower=0> N_size; //
  array[N_size] int<lower=0> year2014;
  array[N_size] int<lower=0> year2015;
  array[N_size] int<lower=0> year2016;
  array[N_size] int<lower=0> year2017;
  array[N_size] int<lower=0> year2018;
  array[N_size] int<lower=0> year2019;
  array[N_size] int<lower=0> year2020;
  array[N_size] int<lower=0> year2021;
  array[N_size] int<lower=0> year2022;
  array[N_size] int<lower=0> year;
  array[N_size] int<lower=0> pyr_on;
  array[N_size] int<lower=0> pyr_off;
  array[N_size,2] int<lower=0> events_linked_on;
  array[N_size,2] int<lower=0> events_linked_off; 
  array[N_size,2] int<lower=0> events_extra; 
  array[N_size] int<lower=0> other_cause_exit; 
  array[N_size] int<lower=0> pyr_died;
  array[N_size] int<lower=0> n_known; 
  array[N_size] int<lower=0> n_known_on;
  array[N_size] int<lower=0> P;
  vector[N_size] age_g2;
  vector[N_size] age_g3;
  vector[N_size] sex_g;
}
transformed data{
  vector<lower=0, upper=1>[N_size] prev_lower_bound;{
    for(i in 1:N_size)
    {
      prev_lower_bound[i] = n_known[i] * 1.0/P[i];
  }
  }
}
//
//
parameters{
  //////////////////////////////////////////////
  vector[2] beta0;
  matrix[12,2] beta;
  vector[3] beta_h;
  vector<lower=0.01, upper=5>[2] sd_ER_re_treat_fy;
  matrix[9,2] EvRate_RE_treat_fy;
  //////////////////////////////////////////////
  vector[12] beta_others;
  //////////////////////////////////////////////
  vector[12] gamma;
  real<lower=0.01, upper=5> sd_PREV_RE_age2_fy;
  real<lower=0.01, upper=5> sd_PREV_RE_age3_fy;
  vector[9] PREV_RE_age2_fy;
  vector[9] PREV_RE_age3_fy;
  //////////////////////////////////////////////
  vector<lower=0.0001,upper=0.1>[N_size] prev_known;
  //////////////////////////////////////////////
}
//
//
//
//
//
//
transformed parameters{
  //Event rates model: off OAT
  matrix[N_size,2] log_lambda_off;
  //Event rates model: on OAT
  matrix[N_size,2] log_lambda_on;
  {
    for(i in 1:N_size){
      for(event in 1:2){
        log_lambda_off[i,event] = beta0[event] + 
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
                          beta_h[2]*age_g3[i]*sex_g[i]*(event-1);
          log_lambda_on[i,event] = log_lambda_off[i,event] + 
                            beta[1,event] +
                            beta_h[3]*sex_g[i]*(event-1) +
                            EvRate_RE_treat_fy[year[i],event];      }
    }
  }
  matrix[N_size,2] lambda_on;
  matrix[N_size,2] lambda_off;
  matrix[N_size,2] mean_events_linked_on;
  matrix[N_size,2] mean_events_linked_off;
  {
    for(i in 1:N_size){
      for(event in 1:2){
        lambda_on[i,event] = exp(log_lambda_on[i,event]);
        mean_events_linked_on[i,event] = lambda_on[i,event] * pyr_on[i];
        lambda_off[i,event] = exp(log_lambda_off[i,event]);
        mean_events_linked_off[i,event] = lambda_off[i,event] * pyr_off[i];
      }
    }
  }
  //////////////////////////////////////////////
  //Event rates model: Other Cause Exit
  vector[N_size] log_lambda_others;
  {
    for(i in 1:N_size){
      log_lambda_others[i] = beta_others[1] + 
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
                              beta_others[12]*year2022[i];
    }
  }
  vector[N_size] lambda_others;
  vector[N_size] mean_other_death;
  {
    for(i in 1:N_size){
      lambda_others[i] = exp(log_lambda_others[i]);
      mean_other_death[i] = lambda_others[i]*pyr_off[i];
    }
  }
  //////////////////////////////////////////////
  //Prevalence model
  vector[N_size] prev;
  vector[N_size] N;
  matrix[N_size,2] mean_events_extra;
  {
  vector[N_size] n_extra;
  vector[N_size] correction_term;
  vector[N_size] pyr_extra;
    for(i in 1:N_size){
      prev[i] = inv_logit(gamma[1] +
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
                        PREV_RE_age3_fy[year[i]]*age_g3[i] 
                        ) + prev_known[i];
      N[i] = P[i] * prev[i];
      n_extra[i] = N[i] - (prev_known[i]*P[i]);
      correction_term[i] = (1/lambda_others[i])*(1 - exp(-lambda_others[i]));
      pyr_extra[i] = pyr_died[i] + (n_extra[i] - pyr_died[i])*correction_term[i]; 
      mean_events_extra[i,1] = lambda_off[i,1] * pyr_extra[i];
      mean_events_extra[i,2] = lambda_off[i,2] * pyr_extra[i];
    }
  }
  //////////////////////////////////////////////
}
//
//
//
//
//
//
model{
  //Likelihood Event Rate Models
  for(i in 1:N_size)
  {
    for(event in 1:2){
      events_linked_off[i,event] ~ poisson(mean_events_linked_off[i,event]);
      events_linked_on[i,event] ~ poisson(mean_events_linked_on[i,event]);
    }
  }
  //Priors for Event Rate Models
  for(i in 1:12)
  {
    for(event in 1:2)
    {
      beta[i,event] ~ normal(0, 100);
    }
  }
  beta0[1] ~ normal(-4.1, 100);
  beta0[2] ~ normal(-4.1, 100);  
  beta_h[1] ~ normal(0, 100);
  beta_h[2] ~ normal(0, 100);
    beta_h[3] ~ normal(0, 100);
  for(y in 1:9){
    for(event in 1:2){
      EvRate_RE_treat_fy[y,event] ~ normal(0,sd_ER_re_treat_fy[event]);
    }
  }
  //////////////////////////////////////////////////
  //Likelihood Other Cause Exit Model
  for(i in 1:N_size){
    other_cause_exit[i] ~ poisson(mean_other_death[i]);
  }
  //Priors for Other Cause Exit Model
  beta_others[1] ~ normal(1, 100);
  for(i in 2:12){
    beta_others[i] ~ normal(0, 100);
  }
  //////////////////////////////////////////////////
  //Likelihood Prevalence Model
  for(i in 1:N_size)
  {
    events_extra[i,1] ~ poisson(mean_events_extra[i,1]); // for deaths
    events_extra[i,2] ~ poisson(mean_events_extra[i,2]); // for hospitalisations
    n_known[i] ~ binomial(P[i],prev_known[i]);
  }
  //Priors for Prevalence Model
  for(i in 1:12){
    gamma ~ normal(0, 100);
  }
  for(y in 1:9){
    PREV_RE_age2_fy[y] ~ normal(0,sd_PREV_RE_age2_fy);
    PREV_RE_age3_fy[y] ~ normal(0,sd_PREV_RE_age3_fy);
  }
  //////////////////////////////////////////////////
}
//
//
//
//
//
//
generated quantities{
////
}

# mpep_model
This project contains example scripts for the implementation of the MPEP model.
It complements the article "Multi-Parameter Estimation of Prevalence (MPEP): A Bayesian modelling approach to estimate the prevalence of opioid dependence", and includes example scripts in Stan and JAGS.

The file "example_data.csv" contains data for estimating the number of people with opioid dependence in Scotland from 2014/15 to 2022/23. The data are stratified by age, sex, and financial year. These data are publicly available (as open data) at the following location:
https://www.opendata.nhs.scot/dataset/estimated-prevalence-of-opioid-dependence-in-scotland/resource/ead97aa5-307d-4d30-a048-3118f2f963fb
No regional stratification is included due to sparse events and restrictions from Public Health Scotland.

The file "model/mpep_model.stan" contains the Stan implementation of the MPEP model.

The file "MPEP_stan.R" includes the R commands to load the data, set initial values, compile the "mpep_model.stan" model, and run the MPEP model.

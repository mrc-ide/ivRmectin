#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
# Create a vector of age categories for the model
#init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 100 #low - 2, moderate - 15, high - 120

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years
time_period <- 730
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

ivm_haz = rep(1, 730) #do this so we can see the new equilibrium
ivm_parms3 <- ivm_fun(IVM_start_times = 180,# time endectocide delivery occurs
                     time_period = time_period,         # time period for the model to run over
                     hazard_profile = ivm_haz[1:730], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                     #hazard_profile = hazzy,
                     ivm_coverage = 0.8, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide

new_mu = 0.1897807
#correctly settles at #28.92708
wh3 <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/endec_mosq_model_check.R",
                                num_int = 3, # number of vector control (IRS and ITN) population groups
                                #het_brackets = 5, # number of heterogeneous biting categories
                                #age = init_age, # the different age classes to be ran within the model
                                init_EIR = init_EIR, # the Entomological Innoculation Rate
                                #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                ttt = ivm_parms3$ttt, # model specific parameter to control timing of endectocide delivery
                                eff_len = ivm_parms3$eff_len, # number of days after receiving endectocide that HR is higher
                                haz = ivm_parms3$haz, # hazard ratio for each off the eff_len number of days
                                ivm_cov_par = ivm_parms3$ivm_cov_par, # proportion of popuulation receiving the endectocide
                                ivm_min_age = ivm_parms3$ivm_min_age, # youngest age group receiving endectocide
                                ivm_max_age = ivm_parms3$ivm_max_age, # oldest age group receiving endectocide
                                IVRM_start = ivm_parms3$IVRM_start,
                                new_mu = new_mu) # model specific parameter to control timing of endectocide delivery
runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}
res3 <- runfun(wh3)
res3$mu_vi
plot(res3$t, res3$mv, ylim = c(0, 50))
res3$mv

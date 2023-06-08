#exploring how feeding behaviour changes with LLIN coverage and does it affect the impact of ivermectin

#using Hannah's odin_model_endectocide.R for this (the full daily hazards model)

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
IVM_start <- c(180, 210, 240)
ivm_cov = 0.8

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}


#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

#IVM on, no nets
mod_1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 1,
                                  #num_int = 2,
                                  #ITN_IRS_on = 100,
                                  #itn_cov = 0.75,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model
res1 <- runfun(mod_1)

mod_2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
                                    het_brackets = 5,
                                    #age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model
res2 <- runfun(mod_2)

plot(res1$t, res1$mv, ylim = c(0, 50))
lines(res2$t, res2$mv, ylim = c(0, 50), col = "red") #ivermectin and nets

#how does av_mosq change over time (avhc parameter = cov[i]*av_mosq)
#av_mosq = rate at which mosquitoes bite each intervention category
#should this be av_human because that is the biting rate on humans in each int category?

#avhc = mean biting rate of mosquitoes in the presence of vector control

plot(res2$t, res2$avhc, ylim = c(0, 1), pch = 1)
#lines(res2$t, res)


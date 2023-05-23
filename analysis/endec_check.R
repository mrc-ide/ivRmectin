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
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

# Load ivermectin hazard
ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
ivm_parms0 <- ivm_fun(IVM_start_times = 10000, #no ivermectin: turning ivermectin on out of bounds of the model run time (3650 days)
                      time_period = time_period,
                      hazard_profile = ivm_haz$IVM_300_3_HS[1:23], #select 300 dosage here
                      #hazard_profile = hazzy,
                      ivm_coverage=ivm_cov,
                      ivm_min_age=5,
                      ivm_max_age = 90)

wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/endec_mosq_model.R",
                                  num_int = 1, #nets only, but no nets actually being rolled out
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms0$ttt,
                                  eff_len = ivm_parms0$eff_len,
                                  haz = ivm_parms0$haz,
                                  ivm_cov_par = ivm_parms0$ivm_cov_par,
                                  ivm_min_age = ivm_parms0$ivm_min_age,
                                  ivm_max_age = ivm_parms0$ivm_max_age,
                                  IVRM_start = ivm_parms0$IVRM_start)

ivm_parms1 <- ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/endec_mosq_model.R",
                                  num_int = 1,
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
res0 <- runfun(wh0)
res1 <- runfun(wh1)

#plotting total number of mosquitoes
cols <- c("grey40", "deeppink2")
plot(res0$t/365, res0$mv, ylim = c(0, 42), main = "Total mosquitoes")
lines(res0$t/365, res1$mv, col = cols[2])
arrows(c(180, 210, 240)/365, -50, c(180, 210, 240)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")

res1_df <- as.data.frame(res1)
res1_IVM_times <- res1_df %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23))
mean_mv_IVM <- mean(res1_IVM_times$mv)
mean_mv_IVM #28.80228

betaa <- unique(res1_df$betaa) #5.489801 #just one betaa, using non-larval model
new_mu = betaa/mean_mv_IVM
new_mu #0.190603

ivm_haz2 = rep(1, 730) #do this so we can see the new equilibrium
ivm_parms3 <- ivm_fun(IVM_start_times = IVM_start,# time endectocide delivery occurs
                     time_period = time_period,         # time period for the model to run over
                     hazard_profile = ivm_haz2[1:23], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                     #hazard_profile = hazzy,
                     ivm_coverage = ivm_cov, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide
new_mu_input = new_mu
#want to settle at mv = 28.80228
wh3 <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/endec_mosq_model_check.R",
                                num_int = 1, # number of vector control (IRS and ITN) population groups
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
                                new_mu = new_mu_input) # model specific parameter to control timing of endectocide delivery

res3 <- runfun(wh3)
res3$mu
plot(res3$t, res3$betaa)
res3$mv[180:730] #off by 1 mosquito, surely this is ok??

plot(res0$t/365, res0$mv, ylim = c(0, 42), main = "Total mosquitoes", xlim = c(0, 3))
lines(res0$t/365, res1$mv, col = cols[2])
abline(h = mean_mv_IVM, lty = "dashed", col = cols[2])
text(x=2, y=mean_mv_IVM+1, 'Mean mv between IVM_start and IVM_end+23')
lines(res0$t/365, res3$mv, col = "cornflowerblue")
arrows(c(180, 210, 240)/365, -50, c(180, 210, 240)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")

#output the prevalence with the full hazard model, 60% coverage and starting IVM on d180

ivm_parms4 <- ivm_fun(IVM_start_times = IVM_start, #no ivermectin: turning ivermectin on out of bounds of the model run time (3650 days)
                      time_period = time_period,
                      hazard_profile = ivm_haz$IVM_300_3_HS[1:23], #select 300 dosage here
                      #hazard_profile = hazzy,
                      ivm_coverage=ivm_cov,
                      ivm_min_age=5,
                      ivm_max_age = 90)

wh4 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 1, #nets only, but no nets actually being rolled out
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms4$ttt,
                                  eff_len = ivm_parms4$eff_len,
                                  haz = ivm_parms4$haz,
                                  ivm_cov_par = ivm_parms4$ivm_cov_par,
                                  ivm_min_age = ivm_parms4$ivm_min_age,
                                  ivm_max_age = ivm_parms4$ivm_max_age,
                                  IVRM_start = ivm_parms4$IVRM_start)
res4 <- runfun(wh4)
#check the mv
plot(res4$t/365, res4$mv, ylim = c(0, 42), main = "Total mosquitoes", xlim = c(0, 3))
#plot the simple model mv onto it
lines(res3$t/365, res3$mv, col = "cornflowerblue")
arrows(c(180, 210, 240)/365, -50, c(180, 210, 240)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")


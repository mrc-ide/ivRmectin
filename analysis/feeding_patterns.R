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
init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

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
                                  age = init_age,
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
                                    age = init_age,
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

#nets only model
#set IVM params for Hannah's mosq model with hazards#
ivm_parms2 <- ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0, # NO IVERMECTIN
  ivm_min_age=5,
  ivm_max_age = 90)

mod_3 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start)

res3 <- runfun(mod_3)


#plot(res1$t, res1$mv, ylim = c(0, 50))
#lines(res2$t, res2$mv, ylim = c(0, 50), col = "red") #ivermectin and nets

#how does av_mosq change over time (avhc parameter = cov[i]*av_mosq)
#av_mosq = rate at which mosquitoes bite each intervention category
#should this be av_human because that is the biting rate on humans in each int category?

#avhc = mean biting rate of mosquitoes in the presence of vector control

plot(res2$t, res2$avhc, ylim = c(0, 1), type = "l")
lines(res2$t, res2$av_mosq[,1], col = "red")
lines(res2$t, res2$av_mosq[,2], col = "blue")


lines(res2$t, res2$av_mosq_sum, col = "turquoise")
lines(res2$t, res2$av_human_sum, col = "maroon")
lines(res2$t, res2$avhc, col = "blue")
lines(res2$t, res2$Q, col = "green")
lines(res2$t, res2$av, col = "red")
lines(res2$t, res2$av_alt_host, col = "purple")


#Hannah saying we might be seeing this behaviour with LLINs because introductions are at a similar time. If we introduce nets earlier..
#and let the model equilibrise at the no-LLIN prevalence (by having a higher starting EIR with the nets), then impact of IVM should be the same

plot(res1$t, res1$slide_prev0to5, ylim = c(0, 1))

time_period22 <- 7000

ivm_parms22 <- ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period22,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

mod_22 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 200,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms22$ttt,
                                    eff_len = ivm_parms22$eff_len,
                                    haz = ivm_parms22$haz,
                                    ivm_cov_par = ivm_parms22$ivm_cov_par,
                                    ivm_min_age = ivm_parms22$ivm_min_age,
                                    ivm_max_age = ivm_parms22$ivm_max_age,
                                    IVRM_start = ivm_parms22$IVRM_start,
                                    ITN_interval = 3*365)
runfun22 <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period22)
  op<- mod$transform_variables(modx)
  return(op)
}
#run Hannah's model
res22 <- runfun22(mod_22)

#plot(res1$t, res1$slide_prev0to5, ylim = c(0, 1), type = "l")
plot(res22$t, res22$slide_prev0to5, col = "blue", type = "l")
lines(res1$t, res1$slide_prev0to5, ylim = c(0, 1), type = "l")

res1$EIR

#are the number of mosquitoes recruited into the ivermectin compartments actually different when nets are on?
plot(res1$t, res1$Sxtot, ylim = c(0, 50), type = "l", ylab = "Susceptible ivm-fed mosquitoes", xlab = "days") # IVM only
lines(res2$t, res2$Sxtot, col = "red") #IVM & LLIN
points(res3$t, res3$Sxtot, col  = "blue") #llin only

#show the prevalence in these two scenarios
plot(res1$t, res1$slide_prev0to5, ylim = c(0, 1), main = "Low residual transmission")
points(res2$t, res2$slide_prev0to5, col  = "red") #ivm and llin
points(res3$t, res3$slide_prev0to5, col  = "blue") #llin only


#impact of ivermectin in an area of high residual malaria transmission
mod_4 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 1,
                                    #num_int = 2,
                                    #ITN_IRS_on = 100,
                                    #itn_cov = 0.75,
                                    #het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start,
                                    bites_Indoors = 0.5,
                                    bites_Bed = 0.5)

#run Hannah's model
res4 <- runfun(mod_4)

#with ivermectin and LLIN
mod_5 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start,
                                    bites_Indoors = 0.5,
                                    bites_Bed = 0.5)

#run Hannah's model
res5 <- runfun(mod_5)

#LLIN only
mod_6 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start,
                                    bites_Indoors = 0.5,
                                    bites_Bed = 0.5)

res6 <- runfun(mod_6)


#high indoor biting
plot(res1$t, res1$slide_prev0to5, ylim = c(0, 1), type = "l", ylab = "Slide prevalence in under 5s", main = "Ivermectin impact under high levels of indoor biting", xlab = "days")
lines(res2$t, res2$slide_prev0to5, col  = "red") #ivm and llin
lines(res3$t, res3$slide_prev0to5, col  = "blue") #llin only
legend(1,1, legend = c("IVM only", "LLIN & IVM", "LLIN only"), col = c("black", "red", "blue"), lty = 1, cex = 0.8)




#low indoor biting
plot(res4$t, res4$slide_prev0to5, col = "green", ylim = c(0, 1), type = "l", ylab = "Slide prevalence in under 5s", main = "Ivermectin impact under low levels of indoor biting", xlab = "days")
lines(res5$t, res5$slide_prev0to5, col = "orange") #ivm and LLIN only, residual malaria
lines(res6$t, res6$slide_prev0to5, col = "purple") #LLIN only
legend(1,1, legend = c("IVM only", "LLIN & IVM", "LLIN only"), col = c("green", "orange", "purple"), lty = 1, cex = 0.8)


plot(res1$t, res1$slide_prev0to5, ylim = c(0, 1), main = "IVM impact depends on residual malaria transmission", type = "l")
lines(res2$t, res2$slide_prev0to5, col  = "red") #ivm and llin
lines(res3$t, res3$slide_prev0to5, col  = "blue") #llin only
lines(res4$t, res4$slide_prev0to5, col = "green") #ivm only, residual malaria
lines(res5$t, res5$slide_prev0to5, col = "orange") #ivm and LLIN only, residual malaria
lines(res6$t, res6$slide_prev0to5, col = "purple") #LLIN only
arrows(c(100, 180, 210, 240), -50, c(100, 180, 210, 240), 0.1, length = 0.1, lwd = 3, col = c ( "cornflowerblue", "goldenrod2", "goldenrod2", "goldenrod2"))
legend(1, 1, legend = c("IVM only. High indoor biting", "IVM & LLIN. High indoor biting",
                         "LLIN only. High indoor biting", "IVM only. Low indoor biting",
                         "IVM & LLIN. Low indoor biting", "LLIN only. Low indoor biting",
                        "LLIN distrib (arrow)", "IVM MDA (arrow)"),
       col =  c("black"  ,"red", "blue", "green", "orange", "purple", "cornflowerblue", "goldenrod2"), lty = 1,  cex = 0.8,)


#look at avhc in model with nets only
plot(res3$t, res3$avhc, ylim = c(0, 1), type = "l")
lines(res3$t, res3$av_alt_host, col = "red")


#different values of Q0
#impact of ivermectin in an area of high residual malaria transmission
mod_7 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    num_int = 1,
                                    #num_int = 2,
                                    #ITN_IRS_on = 100,
                                    #itn_cov = 0.75,
                                    #het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start,
                                   Q0 = 0.7)

#run Hannah's model
res7 <- runfun(mod_7)

#with ivermectin and LLIN
mod_8 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms1$ttt,
                                    eff_len = ivm_parms1$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms1$ivm_cov_par,
                                    ivm_min_age = ivm_parms1$ivm_min_age,
                                    ivm_max_age = ivm_parms1$ivm_max_age,
                                    IVRM_start = ivm_parms1$IVRM_start,
                                    Q0 = 0.7)

#run Hannah's model
res8 <- runfun(mod_8)

#LLIN only
mod_9 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
                                    het_brackets = 5,
                                    age = init_age,
                                    init_EIR = 100,
                                    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                    #admin2 = "Fatick",
                                    ttt = ivm_parms2$ttt,
                                    eff_len = ivm_parms2$eff_len,
                                    haz = ivm_parms1$haz,
                                    ivm_cov_par = ivm_parms2$ivm_cov_par,
                                    ivm_min_age = ivm_parms2$ivm_min_age,
                                    ivm_max_age = ivm_parms2$ivm_max_age,
                                    IVRM_start = ivm_parms2$IVRM_start,
                                    Q0 = 0.7)

res9 <- runfun(mod_9)

plot(res1$t, res1$slide_prev0to5, ylim = c(0, 1), main = "IVM impact and host feeding preference", type = "l")
lines(res2$t, res2$slide_prev0to5, col  = "red") #ivm and llin
lines(res3$t, res3$slide_prev0to5, col  = "blue") #llin only
lines(res7$t, res7$slide_prev0to5, col = "green") #ivm only, residual malaria
lines(res8$t, res8$slide_prev0to5, col = "orange") #ivm and LLIN only, residual malaria
lines(res9$t, res9$slide_prev0to5, col = "purple") #LLIN only
arrows(c(100, 180, 210, 240), -50, c(100, 180, 210, 240), 0.1, length = 0.1, lwd = 3, col = c ( "cornflowerblue", "goldenrod2", "goldenrod2", "goldenrod2"))
legend(1, 1, legend = c("IVM only. High anthropophagy", "IVM & LLIN. High anthropophagy",
                        "LLIN only. High anthropophagy", "IVM only. Moderate anthropophagy",
                        "IVM & LLIN. Moderate anthropophagy", "LLIN only.Moderate anthropophagy",
                        "LLIN distrib (arrow)", "IVM MDA (arrow)"),
       col =  c("black"  ,"red", "blue", "green", "orange", "purple", "cornflowerblue", "goldenrod2"), lty = 1,  cex = 0.8,)

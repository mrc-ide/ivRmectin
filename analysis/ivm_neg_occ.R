#studying the relationship between LLIN usage and IVM killing effect

# x/av_mosq = mu_h(ivm_cov)

#or

#x/avhc = mu_h(ivm_cov)

#avhc = cov*av_mosq

#let's see how av_mosq changes with LLIN coverage
#exploring how feeding behaviour changes with LLIN coverage and does it affect the impact of ivermectin


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

ggplot(ivm_haz, aes(x =))
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

#going to run through a range of ivermectin coverages and track negative mosquito occurrence.####
#not looping because output gets messy
mod_cov_0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                        #num_int = 1,
                                        num_int = 1,
                                        #ITN_IRS_on = 100,
                                        #itn_cov = 0.85,
                                        #het_brackets = 5,
                                        age = init_age,
                                        init_EIR = 100,
                                        #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                        #admin2 = "Fatick",
                                        ttt = ivm_parms1$ttt,
                                        eff_len = ivm_parms1$eff_len,
                                        haz = ivm_parms1$haz,
                                        ivm_cov_par = 0,
                                        ivm_min_age = ivm_parms1$ivm_min_age,
                                        ivm_max_age = ivm_parms1$ivm_max_age,
                                        IVRM_start = ivm_parms1$IVRM_start)


res_cov_0 <- runfun(mod_cov_0)
res_cov_0_df <- as.data.frame(res_cov_0)
neg_cov_0_df <- as.data.frame(which(res_cov_0_df < 0, arr.ind = TRUE))
nrow(neg_cov_0_df) #0 neg occurrences

write.csv(res_cov_0_df, file = "data/res_ivm_cov_0_df.csv", row.names = FALSE)
write.csv(neg_cov_0_df, file = "data/neg_ivm_cov_0_df.csv", row.names = FALSE)

#cov 0.25
mod_cov_0.25 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                           #num_int = 1,
                                           num_int = 1,
                                           #ITN_IRS_on = 100,
                                           #itn_cov = 0.85,
                                           #het_brackets = 5,
                                           age = init_age,
                                           init_EIR = 100,
                                           #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                           #admin2 = "Fatick",
                                           ttt = ivm_parms1$ttt,
                                           eff_len = ivm_parms1$eff_len,
                                           haz = ivm_parms1$haz,
                                           ivm_cov_par = 0.25,
                                           ivm_min_age = ivm_parms1$ivm_min_age,
                                           ivm_max_age = ivm_parms1$ivm_max_age,
                                           IVRM_start = ivm_parms1$IVRM_start)


res_cov_0.25 <- runfun(mod_cov_0.25)
res_cov_0.25_df <- as.data.frame(res_cov_0.25)
neg_cov_0.25_df <- as.data.frame(which(res_cov_0.25_df < 0, arr.ind = TRUE))
nrow(neg_cov_0.25_df) #736 neg occurrences


data %>%
  rowwise() %>%
  filter(any(c_across(starts_with("sam")) > limit))

neg_filt_0.25 <- res_cov_0.25_df[res_cov_0.25_df < 0]
range(neg_filt_0.25)

write.csv(res_cov_0.25_df, file = "data/res_ivm_cov_0.25_df.csv", row.names = FALSE)
write.csv(neg_cov_0.25_df, file = "data/neg_ivm_cov_0.25_df.csv", row.names = FALSE)

res_cov_0.25_df[180, 1475]

#cov 0.5
mod_cov_0.5 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                          #num_int = 1,
                                          num_int = 1,
                                          #ITN_IRS_on = 100,
                                          #itn_cov = 0.85,
                                          #het_brackets = 5,
                                          age = init_age,
                                          init_EIR = 100,
                                          #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                          #admin2 = "Fatick",
                                          ttt = ivm_parms1$ttt,
                                          eff_len = ivm_parms1$eff_len,
                                          haz = ivm_parms1$haz,
                                          ivm_cov_par = 0.5,
                                          ivm_min_age = ivm_parms1$ivm_min_age,
                                          ivm_max_age = ivm_parms1$ivm_max_age,
                                          IVRM_start = ivm_parms1$IVRM_start)


res_cov_0.5 <- runfun(mod_cov_0.5)
res_cov_0.5_df <- as.data.frame(res_cov_0.5)
neg_cov_0.5_df <- as.data.frame(which(res_cov_0.5_df < 0, arr.ind = TRUE))
nrow(neg_cov_0.5_df) #537 neg occurrences

neg_filt_0.5 <- res_cov_0.5_df[res_cov_0.5_df < 0]
range(neg_filt_0.5)

write.csv(res_cov_0.5_df, file = "data/res_ivm_cov_0.5_df.csv", row.names = FALSE)
write.csv(neg_cov_0.5_df, file = "data/neg_ivm_cov_0.5_df.csv", row.names = FALSE)

#cov 0.75
mod_cov_0.75 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                           #num_int = 1,
                                           num_int = 1,
                                           #ITN_IRS_on = 100,
                                           #itn_cov = 0.85,
                                           #het_brackets = 5,
                                           age = init_age,
                                           init_EIR = 100,
                                           #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                           #admin2 = "Fatick",
                                           ttt = ivm_parms1$ttt,
                                           eff_len = ivm_parms1$eff_len,
                                           haz = ivm_parms1$haz,
                                           ivm_cov_par = 0.75,
                                           ivm_min_age = ivm_parms1$ivm_min_age,
                                           ivm_max_age = ivm_parms1$ivm_max_age,
                                           IVRM_start = ivm_parms1$IVRM_start)


res_cov_0.75 <- runfun(mod_cov_0.75)
res_cov_0.75_df <- as.data.frame(res_cov_0.75)
neg_cov_0.75_df <- as.data.frame(which(res_cov_0.75_df < 0, arr.ind = TRUE))
nrow(neg_cov_0.75_df) #643 neg occurrences

neg_filt_0.75 <- res_cov_0.75_df[res_cov_0.75_df < 0]
range(neg_filt_0.75)



write.csv(res_cov_0.75_df, file = "data/res_ivm_cov_0.75_df.csv", row.names = FALSE)
write.csv(neg_cov_0.75_df, file = "data/neg_ivm_cov_0.75_df.csv", row.names = FALSE)

#cov 1
mod_cov_1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                        #num_int = 1,
                                        num_int = 1,
                                        #ITN_IRS_on = 100,
                                        #itn_cov = 0.85,
                                        #het_brackets = 5,
                                        age = init_age,
                                        init_EIR = 100,
                                        #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                        #admin2 = "Fatick",
                                        ttt = ivm_parms1$ttt,
                                        eff_len = ivm_parms1$eff_len,
                                        haz = ivm_parms1$haz,
                                        ivm_cov_par = 1,
                                        ivm_min_age = ivm_parms1$ivm_min_age,
                                        ivm_max_age = ivm_parms1$ivm_max_age,
                                        IVRM_start = ivm_parms1$IVRM_start)


res_cov_1 <- runfun(mod_cov_1)
res_cov_1_df <- as.data.frame(res_cov_1)
neg_cov_1_df <- as.data.frame(which(res_cov_1_df < 0, arr.ind = TRUE))
nrow(neg_cov_1_df) #655 neg occurrences

neg_filt_1 <- res_cov_1_df[res_cov_1_df < 0]
range(neg_filt_1)

x <- filter_all(res_cov_1_df, any_vars(. <0 ))


write.csv(res_cov_1_df, file = "data/res_ivm_cov_1_df.csv", row.names = FALSE)
write.csv(neg_cov_1_df, file = "data/neg_ivm_cov_1_df.csv", row.names = FALSE)

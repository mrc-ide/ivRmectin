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

mu_h_in = 0.257 #this is for 80% IVM coverage

#IVM on, no nets
mod_1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/endec_mosq_model_check.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.75,
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
                                    mu_h = 0)

#run Hannah's model
res1 <- runfun(mod_1)

out <- tail(res1$av_mosq)

plot(res1$t, res1$av_mosq)

mod_2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 0.85,
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
res2 <- runfun(mod_2)

out2 <- tail(res2$av_mosq)

#how much does av_mosq change under different LLIN coverage?

create_ITN_cov_loop <- function(itn_cov_in){
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = 100,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start
  )
  return(output)
}

itn_cov_vector <- seq(0, 1, 0.25)

out_lapply_list <- lapply(itn_cov_vector, create_ITN_cov_loop) #loop through all parameter values

res_out_list <- lapply(out_lapply_list, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot

out_df <- do.call(rbind,
                  sapply(1:length(itn_cov_vector), function(x){
                    as.data.frame(res_out_list[[x]]) %>%
                      select(t, mv, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot) %>%
                      mutate(ref = x)
                  }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(out_df, file = "data/out_df.csv", row.names = FALSE)

#see plot 1 and 2 in llin_ivm_plots.R


#get different combinations of itn cov and llin cov parameters
create_ivm_itn_cov_loop <- function(itn_ivm_param){
  itn_cov_in <- itn_ivm_param[1]
  ivm_cov_in <-itn_ivm_param[2]
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = 100,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start
  )
  return(output)
}

ivm_cov_vector <- seq(0, 1, 0.25)

#make an empty list
param_list <- list()

#use expand grid to get all combinations of the parameter values
param_df <- expand.grid(itn_cov = itn_cov_vector,
                        ivm_cov = ivm_cov_vector)

for(i in seq_len(nrow(param_df))){
  param_list[[i]] <- as.numeric(param_df[i,]) #convert to numeric so it's in the right form for runfun
}

out_2_list <- lapply(param_list, create_ivm_itn_cov_loop) #putting param list into the function to generate parameter set

#run it
res_out_2_list <- lapply(out_2_list, runfun)
#UP TO HERE



#go through and save key parameters


out_df_2 <- do.call(rbind,
                  sapply(1:(length(itn_cov_vector)*length(ivm_cov_vector)), function(x){
                    as.data.frame(res_out_2_list[[x]]) %>%
                      select(t, mv, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot) %>%
                      mutate(ref = x)
                  }, simplify = F))

write.csv(out_df_2, file = "data/out_df_2.csv", row.names = FALSE)




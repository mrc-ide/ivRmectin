#run the model between different values of mu_h. Which value of mu_h gives the best fit to Hannah's model?

require(tidyverse)


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

#set other params#
wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
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
res1 <- runfun(wh1)

#save this output so don't need to run this again and again.
res1_out_df <- as.data.frame(res1)

write.csv(res1_out_df, file = "analysis/exploring_interactions/model_output/endec_model_mv.csv", row.names = FALSE)

#now, trying to find a value of mu_h that best fits Hannah's model

ivm_haz2 = rep(1, 730) #not using hazards, so just set hazards to 1
ivm_parms3 <- ivm_fun(IVM_start_times = IVM_start,# time endectocide delivery occurs
                      #IVM_start = 180,
                      time_period = time_period,         # time period for the model to run over
                      hazard_profile = ivm_haz2[1:23], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                      #hazard_profile = ivm_haz2[1:730],
                      ivm_coverage = ivm_cov, # proportion of population receiving the endectocide
                      ivm_min_age = 5, # youngest age group receiving endectocide
                      ivm_max_age = 90) # oldest age group receiving endectocide

create_mu_h_loop <- function(mu_h_in) {
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide_mu_h.R",
    num_int = 1,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
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
    mu_h = mu_h_in
  )
  return(output)
}

mu_h_vector <- seq(0, 1, 0.001)

#generate the parameter set for the sensitivity analysis
out_lapply_list <- lapply(mu_h_vector, create_mu_h_loop)

#run individual parameters
res_out <- runfun(out_lapply_list[[1]])

#outputting with lapply
res_out_list <- lapply(out_lapply_list, runfun)


res_out_list[[2]]$mu_h #can see outputs for the different mu_h inputs

#read in the model output from file with HRs

res1_out_df_in <- read.csv("analysis/exploring_interactions/model_output/endec_model_mv.csv", header = TRUE)

#select the mvs of Hannah's during IVM distribution time.
#IVM distrib d180, 210 and 240
#so by 240+23, 263, IVM not effective
#so select mv values between 180 and 263.

res1_ivm_distrib <- res1_out_df_in %>%
  filter(between(t, IVM_start[1], ivm_start[3]))


#go through res_out_list (the models with different values of mu_h and extract t, mv and mu_h between t = 180 and t=263)
out_df <- do.call(rbind,
                  sapply(1:length(mu_h_vector), function(x){
                    as.data.frame(res_out_list[[x]]) %>% #go through all the values in the list and save
                      filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
                      select(t, mv, mu_h) %>%
                      mutate(ref=x)
                  }, simplify = F))

write.csv(out_df, file = "data/out_df.csv", row.names = FALSE) #these are the model outputs

out_df_in <- read.csv("data/out_df.csv", header = TRUE)

out_list <- split(out_df_in, f = out_df$ref)

error <- numeric()
for (i in 1:length(mu_h_vector)){
  error <- c(error, sum((res1_ivm_distrib$mv - out_list[[i]]$mv)^2))
}
error
index <- which.min(error) #index the smallest value
mu_h_vector[index] #0.257

mu_val <- 0.257 # be careful with this - might be better to use what is in mu_h_vector[index]http://127.0.0.1:10085/graphics/plot_zoom_png?width=1200&height=900
#mu_val <- 0.43       #with nets

#put this value back into the model
wh4 <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide_mu_h.R",
                                 num_int = 1,
                                 #num_int = 2,# number of vector control (IRS and ITN) population groups
                                 #ITN_IRS_on = 100,
                                 #itn_cov = 0.75,
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
                                 mu_h = mu_val) # model specific parameter to control timing of endectocide delivery

res4 <- runfun(wh4)
res4_out_df <- as.data.frame(res4)
#plots to check output

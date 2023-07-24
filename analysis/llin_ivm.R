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
#library(tidyverse)
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
  modx <- mod$run(t = 1:time_period, atol = 1e-8, rtol = 1e-8)
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
mod_1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 1,
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

df_res1_Sx_F1 <- as.data.frame(res1$Sx_F1)
neg_Sx_F1 <- as.data.frame(which(df_res1_Sx_F1 < 0, arr.ind = TRUE))
min(neg_Sx_F1$row)
df_res1_Sx_F1[180:188,]


out <- tail(res1$av_mosq)

plot(res1$t, res1$av_mosq)

#no nets
mod_1b <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
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
                                     ivm_cov_par = ivm_parms1$ivm_cov_par,
                                     ivm_min_age = ivm_parms1$ivm_min_age,
                                     ivm_max_age = ivm_parms1$ivm_max_age,
                                     IVRM_start = ivm_parms1$IVRM_start)
res_1b <- runfun(mod_1b)
res_1b_df <- as.data.frame(res_1b)
write.csv(res_1b_df, file = "data/res_1b_df.csv", row.names = FALSE)
neg_1b_df <- as.data.frame(which(res_1b_df < 0, arr.ind = TRUE))
nrow(neg_1b_df) #660 negative occurrences
write.csv(neg_1b_df, file = "data/neg_1b_df.csv", row.names = FALSE)
range(res_1b$Sx_F1)
#head(res_1b$Sx_F1)
df_1b_Sx_F1 <- as.data.frame(res_1b$Sx_F1)
neg_Sx_F1 <- as.data.frame(which(df_1b_Sx_F1 <0, arr.ind = TRUE))
min(neg_Sx_F1$row) #in 182
df_1b_Sx_F1[187,]


#trying with a different coverage
mod_1c <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
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


res_1c <- runfun(mod_1c)
res_1c_df <- as.data.frame(res_1c)
neg_1c_df <- as.data.frame(which(res_1c_df < 0, arr.ind = TRUE))
nrow(neg_1c_df) #537 neg occurrences. So it is coverage dependent.

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

write.csv(res_cov_1_df, file = "data/res_ivm_cov_1_df.csv", row.names = FALSE)
write.csv(neg_cov_1_df, file = "data/neg_ivm_cov_1_df.csv", row.names = FALSE)

############
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

range(res2$Sxtot) # 0.000000 4.236532
range(res2$Extot) #goes neg #-3.016024e-14  1.441958e-01
range(res2$Ixtot) #goes neg #-7.533520e-15  3.928892e-03

range(res2$Sx_F1) #-0.0002996689  1.1371498590

df_Sx_F1 <- as.data.frame(res2$Sx_F1)
range(df_Sx_F1$V1) #-4.940656e-324   5.503673e-01
#what row has this value -4.940656e-324

neg_Sx_F1 <- as.data.frame(which(df_Sx_F1 <0, arr.ind = TRUE))
min(neg_Sx_F1$row)
df_Sx_F1[182,] #goes neg in box 3

df_Sx_F2 <- as.data.frame(res2$Sx_F2)
neg_Sx_F2 <- as.data.frame(which(df_Sx_F2 < 0, arr.ind = TRUE))
min(neg_Sx_F2$row) #row 182
df_Sx_F2[182,] #goes neg in box 3

df_Ex_F1 <- as.data.frame(res2$Ex_F1)
df_Ex_F2 <- as.data.frame(res2$Ex_F2)

neg_Ex_F1 <- as.data.frame(which(df_Ex_F1 < 0, arr.ind = TRUE))
min(neg_Ex_F1$row) #row 180
df_Ex_F1[180,] #first neg on d3
neg_Ex_F2 <- as.data.frame(which(df_Ex_F2 < 0, arr.ind = TRUE))
min(neg_Ex_F2$row) #182
df_Ex_F2[182,] #goes neg in box 21

df_Ix_F1 <- as.data.frame(res2$Ix_F1)
df_Ix_F2 <- as.data.frame(res2$Ix_F2)


neg_Ix_F1 <- as.data.frame(which(df_Ix_F1 < 0 , arr.ind = TRUE))
min(neg_Ix_F1$row) #row 276
df_Ix_F1[276,] #box 3
neg_Ix_F2 <- as.data.frame(which(df_Ix_F2 < 0, arr.ind = TRUE))
min(neg_Ix_F2$row) #row 182
df_Ix_F2[182,] #box 3

range(res2$Sx_F2) #-8.232099e-05  4.282043e-01

range(res2$Ex_F1) # -1.279913e-05  1.638044e-02
range(res2$Ex_F2) #-7.631870e-06  9.174634e-03

#just checking neg with Anna's original model

mod_anna <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide_OG.R",
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

res_anna <- runfun(mod_anna)

range(res_anna$Sxtot) # 0.000000 4.236532
range(res_anna$Extot) #-3.016024e-14  1.441958e-01
range(res_anna$Ixtot) #-7.533520e-15  3.928892e-03

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

res_out_2_atol <- lapply(out_2_list, runfun_atol)


#go through and save key parameters


out_df_2 <- do.call(rbind,
                  sapply(1:(length(itn_cov_vector)*length(ivm_cov_vector)), function(x){
                    as.data.frame(res_out_2_list[[x]]) %>%
                     # select(t, mv, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot, FOIv, lag_FOIv) %>%
                      mutate(ref = x)
                  }, simplify = F))

write.csv(out_df_2, file = "data/out_df_2.csv", row.names = FALSE)
require(tidyverse)
out_df_2_atol <- do.call(rbind,
                         sapply(1:(length(itn_cov_vector)*length(ivm_cov_vector)), function(x){
                           as.data.frame(res_out_2_atol[[x]]) %>%
                             #select(t, mv, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot, FOIv, lag_FOIv) %>%
                             mutate(ref = x)
                         }, simplify = F))

write.csv(out_df_2_atol, file = "data/out_df_2_atol.csv", row.names = FALSE)
range(out_df_2$Ex_F1.17)
range(out_df_2_atol$Ex_F1.17)
x <- out_df_2$Ex_F1.17
y <- out_df_2_atol$Ex_F1.17

range(out_df_2$Extot)
range(out_df_2_atol$Extot)

range(out_df_2$Ixtot)
range(out_df_2_atol$Ixtot)

out_df_2$Sx_F1.1[180:250]
table(x == y)

range(out_df_2$Sx_F1.1)
range(out_df_2_atol$Sx_F1.1)

ggplot(out_df_2, aes(x = t, y = Sx_F1.1))+
  geom_point()+
  geom_point(data = out_df_2_atol, aes(x = t, y = Sx_F1.1), col = "red")

#repeating out_df_2, but with low bites_Bed and bites_Indoors (LLINs less impactful here)
create_ivm_itn_cov_loop_evening_vec <- function(itn_ivm_param){
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
    IVRM_start = ivm_parms1$IVRM_start,
    bites_Indoors = 0.5,
    bites_Bed = 0.5
  )
  return(output)
}
itn_cov_vector <- seq(0, 1, 0.25)

ivm_cov_vector <- seq(0, 1, 0.25)

#make an empty list
param_list <- list()

#use expand grid to get all combinations of the parameter values
param_df <- expand.grid(itn_cov = itn_cov_vector,
                        ivm_cov = ivm_cov_vector)

for(i in seq_len(nrow(param_df))){
  param_list[[i]] <- as.numeric(param_df[i,]) #convert to numeric so it's in the right form for runfun
}

out_3_list <- lapply(param_list, create_ivm_itn_cov_loop_evening_vec) #putting param list into the function to generate parameter set

#run it
res_out_3_list <- lapply(out_3_list, runfun)
#UP TO HERE



#go through and save key parameters

require(tidyverse)
out_df_3 <- do.call(rbind,
                    sapply(1:(length(itn_cov_vector)*length(ivm_cov_vector)), function(x){
                      as.data.frame(res_out_3_list[[x]]) %>%
                        select(t, mv, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot, FOIv, lag_FOIv) %>%
                        mutate(ref = x)
                    }, simplify = F))

write.csv(out_df_3, file = "data/out_df_3.csv", row.names = FALSE)

#llin type/resistance and ivermectin uptake#

#for a given net coverage (e.g. 60%, and different ivermectin and r_ITN0 and d_TIN0 values, what is ivermectin uptake?)
#s_ITN0 is estimated from d_ITN0 and r_TIN0

create_ivm_itn_type_vec <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  r_ITN0_in <- itn_type_ivm_param[2]
  #ivm_cov_in <-itn_type_ivm_param[3]
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = 100,
    itn_cov = 0.6,
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
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in
  )
  return(output)
}


runfun_sens <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

d_ITN0_vector <- c(0.25, 0.5, 0.75)
r_ITN0_vector <- 1-d_ITN0_vector
#ivm_cov_vector <- seq(0, 1, 0.25)

#make an empty list
param_list <- list()

#use expand grid to get all combinations of the parameter values
#param_df <- expand.grid(d_ITN0 = d_ITN0_vector,
#                        r_ITN0 = r_ITN0_vector) #,
                        #ivm_cov = ivm_cov_vector)

#can't use expand grid here because d_ITN0 and r_ITN0 are correlated!
param_df <- data.frame(d_ITN0 = d_ITN0_vector, r_ITN0 = r_ITN0_vector)


for(i in seq_len(nrow(param_df))){
  param_list[[i]] <- as.numeric(param_df[i,]) #convert to numeric so it's in the right form for runfun
}

out_4_list <- lapply(param_list, create_ivm_itn_type_vec) #putting param list into the function to generate parameter set

#run it
res_out_4_list <- lapply(out_4_list, runfun_sens)
res_out_4_list[[1]]$d
out_df_4<- do.call(rbind,
                    sapply(1:(length(d_ITN0_vector)), function(x){
                      as.data.frame(res_out_4_list[[x]]) %>%
                        select(t, mv, Q, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot, FOIv, lag_FOIv,
                               d_ITN0, d_ITN, r_ITN0, r_ITN, s_ITN, ITN_decay, itn_loss) %>%
                        mutate(ref = x)
                    }, simplify = F))

write.csv(out_df_4, file = "data/out_df_4.csv", row.names = FALSE)



#end of llin type/resistance and ivermectin uptake#

################################
#now checking ivm coverages (non-zero) with itn cov of 1 to see if getting negative ivermectin mosquitoes
mod_3 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                    #num_int = 1,
                                    num_int = 2,
                                    ITN_IRS_on = 100,
                                    itn_cov = 1,
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

#run Hannah's model
res3 <- runfun(mod_3)
range(res3$Sxtot)
range(res3$Extot)
range(res3$Ixtot)
range(res3$Sv)
range(res3$Ev_F1)
range(res3$Ev_F2)
range(res3$Iv_F1)
range(res3$Iv_F2)

range(res3$Sx_F1) #this is negative
range(res3$Sx_F2) #negative
class(res3$Sx_F1)

class(res3$Sx_F1)
res3_SxF1 <- as.data.frame(res3$Sx_F1)
head(res3_SxF1)

#go from wide to long
res3_SxF1_long <- gather(res3_SxF1, day_feed_after_distrib, number_mosq, V1:V23, factor_key = TRUE)
head(res3_SxF1_long)


res3_SxF1_long <- res3_SxF1_long %>%
  mutate(t = rep(c(seq(1:730)), times = 23),
         is_neg = case_when(number_mosq < 0 ~ "neg",
                            TRUE ~ "pos"))

#plot neg or pos number of mosq, facetted by the fays feeding after distrib
ggplot(res3_SxF1_long, aes(x = t, y = number_mosq,col = as.factor(is_neg)))+
  geom_point()+
  facet_wrap(~day_feed_after_distrib)

neg <- res3_SxF1_long %>%
  filter(is_neg == "neg")
range(neg$t) #187 to 725
unique(neg$day_feed_after_distrib)

head(res3$mu_vi)

res_muvi_long <- as.data.frame(res3$mu_vi) %>%
  gather(day_feed_after_distrib, death_rate, V1:V23, factor_key = TRUE)%>%
  mutate(t = rep(c(seq(1:time_period)), times = ivm_parms1$eff_len))
head(res_muvi_long)

ggplot(res_muvi_long, aes(x = t, y = death_rate, col = as.factor(day_feed_after_distrib)))+
  geom_point()

res3 %>%
  as.data.frame() %>%
  ggplot(aes(x = t, y= mu))+
  geom_point()

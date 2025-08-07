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
library(fBasics)
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

#ggplot(ivm_haz, aes(x =))
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

runfun2 <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, atol = 1e-13, rtol = 1e-14)
  op<- mod$transform_variables(modx)
  return(op)
}#table(res1$Ix_F1<0)[1][[1]]

runfun2a <- function(mod_name, bray1 = 1e-6, bray2 = 1e-6){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, atol = bray1, rtol = bray2)
  op<- mod$transform_variables(modx)
  return(c(16790-table(op$Sx_F1<0)[1][[1]],
           167900-table(op$Ex_F1<0)[1][[1]],
           16790-table(op$Ix_F1<0)[1][[1]],
           min(op$Sx_F1)*Heaviside(-min(op$Sx_F1), a = 0),
           min(op$Ex_F1)*Heaviside(-min(op$Ex_F1), a = 0),
           min(op$Ix_F1)*Heaviside(-min(op$Ix_F1), a = 0)))
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

wh1 <- ivRmectin::create_r_model(
  odin_model_path = "inst/extdata/odin_model_endectocide.R",
  #num_int = 1,
  num_int = 2, # number of vector control (IRS and ITN) population groups
  ITN_IRS_on = 100,
  itn_cov = 0.9,
  #het_brackets = 5, # number of heterogeneous biting categories
  #age = init_age, # the different age classes to be ran within the model
  init_EIR = init_EIR, # the Entomological Innoculation Rate
  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
  #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
  ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
  eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
  haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
  ivm_cov_par = 0.7, # proportion of popuulation receiving the endectocide
  ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
  ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
  IVRM_start = ivm_parms1$IVRM_start
)

res1 <- runfun(wh1)
res2 <- runfun2(wh1)

table(res1$Sx_F1[,17]<0)
table(res2$Sx_F1[,17]<0)
table(res1$Sx_F1[,5]<0)
table(res2$Sx_F1[,5]<0)
table(res1$Sx_F1[,5] == res2$Sx_F1[,5])

table(res1$Sx_F1<0)
table(res2$Sx_F1<0)
#
table(res1$Ex_F1<0)
table(res2$Ex_F1<0)
#
table(res1$Ix_F1<0)
table(res2$Ix_F1<0)

table(res1$Sxtot<0)
table(res2$Sxtot<0)
table(res1$Sxtot == res2$Sxtot)

table(res1$Sv<0)
table(res2$Sv<0)
table(res1$Sv == res2$Sv)
#

storeP <- 10**(seq(-4,-13,-1))
store1_1 <- storeP
store1_2 <- storeP
store1_3 <- storeP
store1_4 <- storeP
store1_5 <- storeP
store1_6 <- storeP
store2_1 <- storeP
store2_2 <- storeP
store2_3 <- storeP
store2_4 <- storeP
store2_5 <- storeP
store2_6 <- storeP
for(i in 1:length(storeP)){
  vb <- runfun2a(wh1, bray1 = storeP[i], bray2 = 10**(-8))
  store1_1[i] <- vb[1]
  store1_2[i] <- vb[2]
  store1_3[i] <- vb[3]
  store1_4[i] <- vb[4]
  store1_5[i] <- vb[5]
  store1_6[i] <- vb[6]
  vb2 <- runfun2a(wh1, bray2 = storeP[i])
  store2_1[i] <- vb2[1]
  store2_2[i] <- vb2[2]
  store2_3[i] <- vb2[3]
  store2_4[i] <- vb2[4]
  store2_5[i] <- vb2[5]
  store2_6[i] <- vb2[6]
}

df1 <- data.frame('atol' = log10(storeP),'x1' = store1_1,'x2' = store1_2,
                  'x3' = store1_3,'x4' = store1_4,
                  'x5' = store1_5,'x6' = store1_6)
df2 <- data.frame('rtol' = log10(storeP),'x1' = store2_1,'x2' = store2_2,
                  'x3' = store2_3,'x4' = store2_4,
                  'x5' = store2_5,'x6' = store2_6)
df1m <- reshape2::melt(df1,id.vars = 'precision')
#ggplot(df1m) + geom_line(aes(x=precision, y=value, color = variable)) + theme_classic()
ggplot(df1) + geom_line(aes(x=atol, y=x1)) + theme_classic()
ggplot(df1) + geom_line(aes(x=atol, y=x2)) + theme_classic()
ggplot(df1) + geom_line(aes(x=atol, y=x3)) + theme_classic()
#
ggplot(df1) + geom_line(aes(x=atol, y=log10(abs(x4)))) + theme_classic()
ggplot(df1) + geom_line(aes(x=atol, y=log10(abs(x5)))) + theme_classic()
ggplot(df1) + geom_line(aes(x=atol, y=log10(abs(x6)))) + theme_classic()


ggplot(df2) + geom_line(aes(x=rtol, y=x1)) + theme_classic()
ggplot(df2) + geom_line(aes(x=rtol, y=x2)) + theme_classic()
ggplot(df2) + geom_line(aes(x=rtol, y=x3)) + theme_classic()
#
ggplot(df2) + geom_line(aes(x=rtol, y=log10(abs(x4)))) + theme_classic()
ggplot(df2) + geom_line(aes(x=rtol, y=log10(abs(x5)))) + theme_classic()
ggplot(df2) + geom_line(aes(x=rtol, y=log10(abs(x6)))) + theme_classic()



























#
#
#
#
# #get different combinations of itn cov and llin cov parameters
# create_ivm_itn_cov_loop <- function(itn_ivm_param){
#   itn_cov_in <- itn_ivm_param[1]
#   ivm_cov_in <-itn_ivm_param[2]
#   output <- ivRmectin::create_r_model(
#     odin_model_path = "inst/extdata/odin_model_endectocide.R",
#     #num_int = 1,
#     num_int = 2, # number of vector control (IRS and ITN) population groups
#     ITN_IRS_on = 100,
#     itn_cov = itn_cov_in,
#     #het_brackets = 5, # number of heterogeneous biting categories
#     #age = init_age, # the different age classes to be ran within the model
#     init_EIR = init_EIR, # the Entomological Innoculation Rate
#     #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
#     #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
#     ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
#     eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
#     haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
#     ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
#     ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
#     ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
#     IVRM_start = ivm_parms1$IVRM_start
#   )
#   return(output)
# }
#
# itn_cov_vector <- seq(0, 1, 0.25)
# ivm_cov_vector <- seq(0, 1, 0.25)
#
# #make an empty list
# param_list <- list()
#
# #use expand grid to get all combinations of the parameter values
# param_df <- expand.grid(itn_cov = itn_cov_vector,
#                         ivm_cov = ivm_cov_vector)
#
# for(i in seq_len(nrow(param_df))){
#   param_list[[i]] <- as.numeric(param_df[i,]) #convert to numeric so it's in the right form for runfun
# }
#
# out_2_list <- lapply(param_list, create_ivm_itn_cov_loop) #putting param list into the function to generate parameter set
#
# #run it
# res_out_2_list <- lapply(out_2_list, runfun)
# #UP TO HERE
#
# res_out_2_atol <- lapply(out_2_list, runfun_atol)
#
#
# #go through and save key parameters
#
# require(tidyverse) #loading this in early can mask functions in ivRmectin
# out_df_2 <- do.call(rbind,
#                     sapply(1:(length(itn_cov_vector)*length(ivm_cov_vector)), function(x){
#                       as.data.frame(res_out_2_list[[x]]) %>%
#                         # select(t, mv, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot, FOIv, lag_FOIv) %>%
#                         mutate(ref = x)
#                     }, simplify = F))
#
# write.csv(out_df_2, file = "data/out_df_2.csv", row.names = FALSE)
# require(tidyverse)
# out_df_2_atol <- do.call(rbind,
#                          sapply(1:(length(itn_cov_vector)*length(ivm_cov_vector)), function(x){
#                            as.data.frame(res_out_2_atol[[x]]) %>%
#                              #select(t, mv, avhc, av_mosq_sum, itn_cov, ivm_cov, Sxtot, Extot, Ixtot, mvxtot, FOIv, lag_FOIv) %>%
#                              mutate(ref = x)
#                          }, simplify = F))
#
# write.csv(out_df_2_atol, file = "data/out_df_2_atol.csv", row.names = FALSE)

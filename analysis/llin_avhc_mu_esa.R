#Do LLINs influence the excess mortality rate due to ivermectin#?
#going to use this analysis to see if we need to change the way in which we model excess mortality due to ivm in malariasimulation

#odin_endec_mu_h_model.R has the mu_vi[i] = haz[i]*(mu+mu_h), where haz is always 1

#odin_endec_mu_h_constant_ivm_uptake.R #this model has a constant parameter for the uptake of ivermectin

require(tidyverse)


#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
require(tidyverse)
# Create a vector of age categories for the model
#init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 100 #low - 2, moderate - 15, high - 120

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years
time_period <- 365*10
IVM_begin <- (365*8)+180
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)
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

eff_len <- 23
ivm_on <- IVM_start[1]
ivm_off <- IVM_start[3]+eff_len

ivm_haz2 = rep(1, time_period) #not using hazards, so just set hazards to 1
ivm_parms3 <- ivm_fun(IVM_start_times = IVM_start,# time endectocide delivery occurs
                      #IVM_start = 180,
                      time_period = time_period,         # time period for the model to run over
                      hazard_profile = ivm_haz2[1:23], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                      #hazard_profile = ivm_haz2[1:730],
                      ivm_coverage = ivm_cov, # proportion of population receiving the endectocide
                      ivm_min_age = 5, # youngest age group receiving endectocide
                      ivm_max_age = 90) # oldest age group receiving endectocide

#run Hannah's model, output the avhc and mosquito density over time
wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  #num_int = 1,
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

#run Hannah's model
res0 <- runfun(wh0)
df_0 <- as.data.frame(res0)

#NC model
#ivm only, no nets
wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_endec_mu_h_model.R",
                                  #num_int = 1,
                                  num_int = 1,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms3$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  mu_h = 0.26,
                                  IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model with mu_h of 0.26
res2 <- runfun(wh2)
df_2 <- as.data.frame(res2)

#df_0 : Hannah's model, no nets, ivm
df_0_epi <- df_0 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead)
write.csv(df_0_epi, file = "data/llin_ivm_muh/df_0_epi.csv", row.names = FALSE)

#df2: Nilani's model, no nets, ivm
df_2_epi <- df_2 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead)
write.csv(df_2_epi, file = "data/llin_ivm_muh/df_2_epi.csv", row.names = FALSE)

#get net eff profiles for different net types
#going to look at resistance 10%, 50%, 70%, 90%
#and also Ettie's UTNs: ettie_rn0 = 0.409, ettie_dn0 = 0.059

res_vector <- c(0, 0.1, 0.5, 0.7, 0.9)

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance == 0 | resistance  == 0.1 | resistance == 0.5 | resistance == 0.7 | resistance == 0.9)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

IG2_d_ITN0 <- df_IG2$dn0_med
pyr_only_d_ITN0 <- df_pyr_only$dn0_med
pyr_pbo_d_ITN0 <- df_pyr_pbo$dn0_med


pyr_param_list <- list()
pyr_pbo_param_list <- list()
IG2_param_list <- list()


itn_cov_vector <- seq(0, 0.8, 0.2)

#only use expand grid on the ivm_cov and dn0, because rn0 is correlated with dn0 so we add it in via left_join
IG2_param_df_crit <- expand.grid(dn0_med = IG2_d_ITN0, itn_cov = itn_cov_vector) #606 rows
require(tidyverse)
IG2_param_df <- left_join(IG2_param_df_crit,
                          df_IG2 %>% dplyr::select(dn0_med, rn0_med),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "IG2") %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med)

#repeat for pyrethroid nets
pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vector)
pyr_param_df <- left_join(pyr_param_df_crit,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med)
#repeat for pbo nets
pbo_param_df_crit <- expand.grid(dn0_med = pyr_pbo_d_ITN0, itn_cov = itn_cov_vector)
pbo_param_df <- left_join(pbo_param_df_crit,
                          df_pyr_pbo %>% dplyr::select(dn0_med, rn0_med),
                          by = c("dn0_med")) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med)


#utn parameters
#ettie_rn0 = 0.409, ettie_dn0 = 0.059
ettie_dn0 = 0.059
ettie_rn0 = 0.409

df_utn <- data.frame(ettie_dn0, ettie_rn0, "resistance" ,"UTN")
names(df_utn) <- c("d_ITN0", "r_ITN0", "99" ,"net_type")
utn_param_df <- expand.grid(dn0_med = df_utn$d_ITN0, itn_cov = itn_cov_vector) %>%
  mutate(r_ITN0 = df_utn$r_ITN0) %>%
  rename(d_ITN0 = dn0_med)

df_utn_list <- list()

#conversions to numeric
for (i in seq_len(nrow(utn_param_df))){
  df_utn_list[[i]] <- as.numeric(utn_param_df[i,])
}

for (i in seq_len(nrow(pyr_param_df))){
  pyr_param_list[[i]] <- as.numeric(pyr_param_df[i,])
}

for(i in seq_len(nrow(IG2_param_df))){
  IG2_param_list[[i]] <- as.numeric(IG2_param_df[i,])
}


for (i in seq_len(nrow(pbo_param_df))){
  pyr_pbo_param_list[[i]] <- as.numeric(pbo_param_df[i,])
}


HS_ITN_cov_loop <- function(itn_type_param){
  d_ITN0_in <- itn_type_param[1]
  itn_cov_in <- itn_type_param[2]
  r_ITN0_in <- itn_type_param[3]
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
    ivm_cov_par = ivm_parms1$ivm_cov_par, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in
  )
  return(output)
}



utn_out_list <- lapply(df_utn_list, HS_ITN_cov_loop) #loop through all parameter values

res_utn_out <- lapply(utn_out_list, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
utn_out_df_HS <- do.call(rbind,
                     sapply(1:length(itn_cov_vector), function(x){
                       as.data.frame(res_utn_out[[x]]) %>%
                         select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                d_ITN0, r_ITN0) %>%
                         mutate(ref = x, net_type = "UTN")
                     }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(utn_out_df_HS, file = "data/llin_ivm_muh/utn_out_df_HS.csv", row.names = FALSE)


#pyrethroids

pyr_out_list <- lapply(pyr_param_list, HS_ITN_cov_loop) #loop through all parameter values

res_pyr_out <- lapply(pyr_out_list, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
pyr_out_df_HS <- do.call(rbind,
                      sapply(1:(length(itn_cov_vector)*length(res_vector)), function(x){
                        as.data.frame(res_pyr_out[[x]]) %>%
                          select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                 d_ITN0, r_ITN0) %>%
                          mutate(ref = x, net_type = "pyr_only")
                      }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(pyr_out_df_HS, file = "data/llin_ivm_muh/pyr_out_df_HS.csv", row.names = FALSE)

#IG2s

IG2_out_list <- lapply(IG2_param_list, HS_ITN_cov_loop) #loop through all parameter values

res_IG2_out <- lapply(IG2_out_list, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
IG2_out_df_HS <- do.call(rbind,
                         sapply(1:(length(itn_cov_vector)*length(res_vector)), function(x){
                           as.data.frame(res_IG2_out[[x]]) %>%
                             select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                    d_ITN0, r_ITN0) %>%
                             mutate(ref = x, net_type = "IG2")
                         }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(IG2_out_df_HS, file = "data/llin_ivm_muh/IG2_out_df_HS.csv", row.names = FALSE)

#pyr pbo pbo_param_df
pyr_pbo_out_list <- lapply(pyr_pbo_param_list, HS_ITN_cov_loop) #loop through all parameter values

res_pyr_pbo_out <- lapply(pyr_pbo_out_list, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
pyr_pbo_out_df_HS <- do.call(rbind,
                         sapply(1:(length(itn_cov_vector)*length(res_vector)), function(x){
                           as.data.frame(res_pyr_pbo_out[[x]]) %>%
                             select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                    d_ITN0, r_ITN0) %>%
                             mutate(ref = x, net_type = "pyr_pbo")
                         }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(pyr_pbo_out_df_HS, file = "data/llin_ivm_muh/pyr_pbo_out_df_HS.csv", row.names = FALSE)



avhc_constant <- 0.3066667

NC_ITN_cov_loop <- function(itn_type_param){
  d_ITN0_in <- itn_type_param[1]
  itn_cov_in <- itn_type_param[2]
  r_ITN0_in <- itn_type_param[3]
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_endec_mu_h_constant_ivm_uptake.R",
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
    haz = ivm_parms3$haz, # hazard ratio for each off the eff_len number of days. SET TO 1
    ivm_cov_par = ivm_parms1$ivm_cov_par, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    avhc = avhc_constant,
    mu_h = 0.26, #from the fitting
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in
  )
  return(output)
}

utn_out_list_NC <- lapply(df_utn_list, NC_ITN_cov_loop) #loop through all parameter values

res_utn_out_NC <- lapply(utn_out_list_NC, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
utn_out_df_NC <- do.call(rbind,
                         sapply(1:length(itn_cov_vector), function(x){
                           as.data.frame(res_utn_out_NC[[x]]) %>%
                             select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                    d_ITN0, r_ITN0) %>%
                             mutate(ref = x, net_type = "UTN")
                         }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(utn_out_df_NC, file = "data/llin_ivm_muh/utn_out_df_NC.csv", row.names = FALSE)


#pyrethroids

pyr_out_list_NC <- lapply(pyr_param_list, NC_ITN_cov_loop) #loop through all parameter values

res_pyr_out_NC <- lapply(pyr_out_list_NC, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
pyr_out_df_NC <- do.call(rbind,
                         sapply(1:(length(itn_cov_vector)*length(res_vector)), function(x){
                           as.data.frame(res_pyr_out_NC[[x]]) %>%
                             select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                    d_ITN0, r_ITN0) %>%
                             mutate(ref = x, net_type = "pyr_only")
                         }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(pyr_out_df_NC, file = "data/llin_ivm_muh/pyr_out_df_NC.csv", row.names = FALSE)

#IG2s

IG2_out_list_NC <- lapply(IG2_param_list, NC_ITN_cov_loop) #loop through all parameter values

res_IG2_out_NC <- lapply(IG2_out_list_NC, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
IG2_out_df_NC <- do.call(rbind,
                         sapply(1:(length(itn_cov_vector)*length(res_vector)), function(x){
                           as.data.frame(res_IG2_out_NC[[x]]) %>%
                             select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                    d_ITN0, r_ITN0) %>%
                             mutate(ref = x, net_type = "IG2")
                         }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(IG2_out_df_NC, file = "data/llin_ivm_muh/IG2_out_df_NC.csv", row.names = FALSE)

#pyr pbo pbo_param_df
pyr_pbo_out_list_NC <- lapply(pyr_pbo_param_list, NC_ITN_cov_loop) #loop through all parameter values

res_pyr_pbo_out_NC <- lapply(pyr_pbo_out_list_NC, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
pyr_pbo_out_df_NC <- do.call(rbind,
                             sapply(1:(length(itn_cov_vector)*length(res_vector)), function(x){
                               as.data.frame(res_pyr_pbo_out_NC[[x]]) %>%
                                 select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead,
                                        d_ITN0, r_ITN0) %>%
                                 mutate(ref = x, net_type = "pyr_pbo")
                             }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(pyr_pbo_out_df_NC, file = "data/llin_ivm_muh/pyr_pbo_out_df_NC.csv", row.names = FALSE)

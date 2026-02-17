#model runs for each scenario

#script to compare different model types
devtools::load_all()
require(tidyverse)

#get model in right space for each setting: get EIR for each setting

years <- seq(2013, 2025, by = 1)
distr_years <- c(years[6], years[9], years[13])
#say ITNs are distributed in March
distr_campaign <- c(30*3)/365
distr_campaign_days <- distr_campaign*365

time_period <- 365*13

net_seq <- seq(distr_campaign_days+(365*5), time_period, by = 3*365)
itn_on <- net_seq[1]

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#prevalence survey is on 1st Jan 2024
y_2024 <- 365*12
prev_survey_date <- 1+y_2024 #1st Jan 2024


#start MDA in March
mda_start <- (30*3)
IVM_begin1 <- (365*12)+mda_start

mda_int <- 30 #every 30 days

IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

ivm_IRR<- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/IVER_ONE_dat/sequential_model_irr_table.csv", header = TRUE)

colnames(ivm_IRR) = c("Day", "IVM_400_1", "IVM_800_1", "IVM_300_3", "IVM_600_3")

ivm_IRR <- ivm_IRR %>%
  filter(Day !=0)

# Sourcing the extra functions required to generate the endec parameters
source("R/mda_ivm_functions.R")

#set the IVM parms for the different dosages
#1x400 is above 1 until day 58
#1x800 until day 60
#3x300 until day 60
#3x600 until day 60


ivm_parms1_400 <- ivRmectin::ivm_fun(
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_IRR$IVM_400_1[1:58],
  ivm_coverage=0.7)

ivm_parms1_800 <- ivRmectin::ivm_fun(
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_IRR$IVM_800_1[1:60],
  ivm_coverage=0.7)

ivm_parms1_3_300 <- ivRmectin::ivm_fun(
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_IRR$IVM_300_3[1:60],
  ivm_coverage=0.7)

ivm_parms1_3_600 <- ivRmectin::ivm_fun(
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_IRR$IVM_600_3[1:60],
  ivm_coverage=0.7)


####KWALE####
param_df_kwale <-readRDS("analysis/IVER_ONE/kwale_param.rds")
explore_kwale <- param_df_kwale %>%
  select(-(init_EIR))

#want to hit ~12% u5 prevalence on 1st Jan 2024

init_EIR_vec_kwale <- c(10,12,15)
mod_space_kwale <- explore_kwale_expanded <- merge(
  explore_kwale,
  expand.grid(init_EIR = init_EIR_vec_kwale),
  by = NULL
)

mod_space_kwale_list <- list()

for(i in seq_len(nrow(mod_space_kwale))){
  mod_space_kwale_list[[i]] <- as.numeric(mod_space_kwale[i,])
}

ento_space_mod_kwale <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_all_age.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Kwale", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_300$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_300$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_300$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on
  )
  return(output)
}

my_sim_kwale <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_space_kwale_list, ento_space_mod_kwale)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_space_kwale)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "explore", int = "ITN-only"))}, simplify = F))
  return(out_df_antag)

}

my_sim_kwale_out <- my_sim_kwale()


ggplot(my_sim_kwale_out, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  facet_wrap(vars(ref))+
  geom_vline(xintercept = prev_survey_date)

my_sim_kwale_out %>%
  group_by(ref) %>%
  filter(t == prev_survey_date) %>%
  summarise(prev = slide_prev0to5)

#ref 2 is pretty good, that is an init_EIR of 12
params_kwale <- mod_space_kwale %>%
  filter(init_EIR == init_EIR_vec_kwale[2])

#BUSIA####
param_df_busia <-readRDS("analysis/IVER_ONE/busia_param.rds")
explore_busia <- param_df_busia %>%
  select(-(init_EIR))

init_EIR_vec_busia <- c(30,32,34)
mod_space_busia <- explore_busia_expanded <- merge(
  explore_busia,
  expand.grid(init_EIR = init_EIR_vec_busia),
  by = NULL
)

mod_space_busia_list <- list()

for(i in seq_len(nrow(mod_space_busia))){
  mod_space_busia_list[[i]] <- as.numeric(mod_space_busia[i,])
}

ento_space_mod_busia <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_all_age.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Busia", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_300$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_300$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_300$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on
  )
  return(output)
}

my_sim_busia <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_space_busia_list, ento_space_mod_busia)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_space_busia)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "explore", int = "ITN-only"))}, simplify = F))
  return(out_df_antag)

}

my_sim_busia_out <- my_sim_busia()


ggplot(my_sim_busia_out, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  facet_wrap(vars(ref))+
  geom_vline(xintercept = prev_survey_date)

my_sim_busia_out %>%
  group_by(ref) %>%
  filter(t == prev_survey_date) %>%
  summarise(prev = slide_prev0to5)

#ref 3 EIR close enough to 22%

params_busia <- mod_space_busia %>%
  filter(init_EIR == init_EIR_vec_busia[3])

#HOMA BAY####
param_df_homa_bay <-readRDS("analysis/IVER_ONE/homa_bay_param.rds")
explore_homa_bay <- param_df_homa_bay %>%
  select(-(init_EIR))

init_EIR_vec_homa_bay <- c(1,2,3)
mod_space_homa_bay <- explore_homa_bay_expanded <- merge(
  explore_homa_bay,
  expand.grid(init_EIR = init_EIR_vec_homa_bay),
  by = NULL
)

mod_space_homa_bay_list <- list()

for(i in seq_len(nrow(mod_space_homa_bay))){
  mod_space_homa_bay_list[[i]] <- as.numeric(mod_space_homa_bay[i,])
}

ento_space_mod_homa_bay <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_all_age.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Homa Bay", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_300$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_300$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_300$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on
  )
  return(output)
}

my_sim_homa_bay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_space_homa_bay_list, ento_space_mod_homa_bay)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_space_homa_bay)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "explore", int = "ITN-only"))}, simplify = F))
  return(out_df_antag)

}

my_sim_homa_bay_out <- my_sim_homa_bay()

ggplot(my_sim_homa_bay_out, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  facet_wrap(vars(ref))+
  geom_vline(xintercept = prev_survey_date)

#find close to 1.1% prev u5

my_sim_homa_bay_out %>%
  group_by(ref) %>%
  filter(t == prev_survey_date) %>%
  summarise(prev = slide_prev0to5)

#ref 1 EIR close enough to 1.1% prev but be very careful - this is super low!

params_homa_bay <- mod_space_homa_bay %>%
  filter(init_EIR == init_EIR_vec_homa_bay[1])

#MIGORI####
param_df_migori <-readRDS("analysis/IVER_ONE/migori_param.rds")
explore_migori <- param_df_migori %>%
  select(-(init_EIR))

init_EIR_vec_migori <- c(8,9,10)

mod_space_migori <- explore_migori <- merge(
  explore_migori,
  expand.grid(init_EIR = init_EIR_vec_migori),
  by = NULL
)

mod_space_migori_list <- list()

for(i in seq_len(nrow(mod_space_migori))){
  mod_space_migori_list[[i]] <- as.numeric(mod_space_migori[i,])
}

ento_space_mod_migori <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_all_age.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Migori", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_300$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_300$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_300$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on
  )
  return(output)
}

my_sim_migori <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_space_migori_list, ento_space_mod_migori)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_space_migori)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "explore", int = "ITN-only"))}, simplify = F))
  return(out_df_antag)

}

my_sim_migori_out <- my_sim_migori()

ggplot(my_sim_migori_out, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  facet_wrap(vars(ref))+
  geom_vline(xintercept = prev_survey_date)

#find close to 9% prev u5

my_sim_migori_out %>%
  group_by(ref) %>%
  filter(t == prev_survey_date) %>%
  summarise(prev = slide_prev0to5)

#ref 2 EIR close enough to 9% prev

params_migori <- mod_space_migori %>%
  filter(init_EIR == init_EIR_vec_migori[2])

#SIAYA####

param_df_siaya <-readRDS("analysis/IVER_ONE/siaya_param.rds")
explore_siaya <- param_df_siaya %>%
  select(-(init_EIR))

init_EIR_vec_siaya <- c(35,40,45)

mod_space_siaya <- explore_siaya <- merge(
  explore_siaya,
  expand.grid(init_EIR = init_EIR_vec_siaya),
  by = NULL
)

mod_space_siaya_list <- list()

for(i in seq_len(nrow(mod_space_siaya))){
  mod_space_siaya_list[[i]] <- as.numeric(mod_space_siaya[i,])
}

ento_space_mod_siaya <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_all_age.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Siaya", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_300$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_300$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_300$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on
  )
  return(output)
}

my_sim_siaya <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_space_siaya_list, ento_space_mod_siaya)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_space_siaya)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "explore", int = "ITN-only"))}, simplify = F))
  return(out_df_antag)

}

my_sim_siaya_out <- my_sim_siaya()

ggplot(my_sim_siaya_out, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  facet_wrap(vars(ref))+
  geom_vline(xintercept = prev_survey_date)

#find close to 25% prev u5

my_sim_siaya_out %>%
  group_by(ref) %>%
  filter(t == prev_survey_date) %>%
  summarise(prev = slide_prev0to5)

#ref 2 EIR close enough to 25% prev

params_siaya <- mod_space_siaya %>%
  filter(init_EIR == init_EIR_vec_siaya[2])

#save the outputs
saveRDS(params_kwale, file = "analysis/IVER_ONE/param_input_odin_ivm_model/params_kwale.rds")
saveRDS(params_busia, file = "analysis/IVER_ONE/param_input_odin_ivm_model/params_busia.rds")
saveRDS(params_homa_bay, file = "analysis/IVER_ONE/param_input_odin_ivm_model/params_homa_bay.rds")
saveRDS(params_migori, file = "analysis/IVER_ONE/param_input_odin_ivm_model/params_migori.rds")
saveRDS(params_siaya, file = "analysis/IVER_ONE/param_input_odin_ivm_model/params_siaya.rds")





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


#read in the different param sets
mod_df_kwale <- readRDS("analysis/IVER_ONE/param_input_odin_ivm_model/params_kwale.rds")
mod_df_busia <- readRDS("analysis/IVER_ONE/param_input_odin_ivm_model/params_busia.rds")
mod_df_homa_bay <- readRDS("analysis/IVER_ONE/param_input_odin_ivm_model/params_homa_bay.rds")
mod_df_migori <- readRDS("analysis/IVER_ONE/param_input_odin_ivm_model/params_migori.rds")
mod_df_siaya <- readRDS("analysis/IVER_ONE/param_input_odin_ivm_model/params_siaya.rds")

mod_kwale_list <- list()
for(i in seq_len(nrow(mod_df_kwale))){
  mod_kwale_list[[i]] <- as.numeric(mod_df_kwale[i,])
}

mod_busia_list <- list()
for(i in seq_len(nrow(mod_df_busia))){
  mod_busia_list[[i]] <- as.numeric(mod_df_busia[i,])
}

mod_homa_bay_list <- list()
for(i in seq_len(nrow(mod_df_homa_bay))){
  mod_homa_bay_list[[i]] <- as.numeric(mod_df_homa_bay[i,])
}

mod_migori_list <- list()
for(i in seq_len(nrow(mod_df_migori))){
  mod_migori_list[[i]] <- as.numeric(mod_df_migori[i,])
}

mod_siaya_list <- list()
for(i in seq_len(nrow(mod_df_siaya))){
  mod_siaya_list[[i]] <- as.numeric(mod_df_siaya[i,])
}



#KWALE####
#1x400
ento_mod_kwale_1_400 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
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

my_sim_kwale_1_400 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_kwale_list, ento_mod_kwale_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_kwale)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "1_400"))}, simplify = F))
  return(out_df_antag)

}

my_sim_kwale_out_1_400 <- my_sim_kwale_1_400()

saveRDS(my_sim_kwale_out_1_400, file = "analysis/IVER_ONE/output/out_kwale_1_400.rds")

ggplot(my_sim_kwale_out_1_400, aes(x = t, y = EIRout))+
  geom_line()

#1x800
ento_mod_kwale_1_800 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
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

my_sim_kwale_1_800 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_kwale_list, ento_mod_kwale_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_kwale)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "1_800"))}, simplify = F))
  return(out_df_antag)

}

my_sim_kwale_out_1_800 <- my_sim_kwale_1_800()

saveRDS(my_sim_kwale_out_1_800, file = "analysis/IVER_ONE/output/out_kwale_1_800.rds")

ggplot(my_sim_kwale_out_1_800, aes(x = t, y = EIRout))+
  geom_line()

#3x300
ento_mod_kwale_3_300 <- function(itn_type_ivm_param){
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
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

my_sim_kwale_3_300 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_kwale_list, ento_mod_kwale_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_kwale)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "3_300"))}, simplify = F))
  return(out_df_antag)

}

my_sim_kwale_out_3_300 <- my_sim_kwale_3_300()

saveRDS(my_sim_kwale_out_3_300, file = "analysis/IVER_ONE/output/out_kwale_3_300.rds")

ggplot(my_sim_kwale_out_3_300, aes(x = t, y = EIRout))+
  geom_line()

#3x600
ento_mod_kwale_3_600 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
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

my_sim_kwale_3_600 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_kwale_list, ento_mod_kwale_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_kwale)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "3_600"))}, simplify = F))
  return(out_df_antag)

}

my_sim_kwale_out_3_600 <- my_sim_kwale_3_600()

saveRDS(my_sim_kwale_out_3_600, file = "analysis/IVER_ONE/output/out_kwale_3_600.rds")

ggplot(my_sim_kwale_out_3_300, aes(x = t, y = EIRout))+
  geom_line()+
  geom_line(data = my_sim_kwale_out_3_600, aes(x = t, y = EIRout), col = "red")


#BUSIA#####

#1x400
ento_mod_busia_1_400 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
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

my_sim_busia_1_400 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_busia_list, ento_mod_busia_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_busia)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "1_400"))}, simplify = F))
  return(out_df_antag)

}

my_sim_busia_out_1_400 <- my_sim_busia_1_400()

saveRDS(my_sim_busia_out_1_400, file = "analysis/IVER_ONE/output/out_busia_1_400.rds")

#1x800
ento_mod_busia_1_800 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
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

my_sim_busia_1_800 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_busia_list, ento_mod_busia_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_busia)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "1_800"))}, simplify = F))
  return(out_df_antag)

}

my_sim_busia_out_1_800 <- my_sim_busia_1_800()

saveRDS(my_sim_busia_out_1_800, file = "analysis/IVER_ONE/output/out_busia_1_800.rds")

#3x300
ento_mod_busia_3_300 <- function(itn_type_ivm_param){
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
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

my_sim_busia_3_300 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_busia_list, ento_mod_busia_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_busia)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "3_300"))}, simplify = F))
  return(out_df_antag)

}

my_sim_busia_out_3_300 <- my_sim_busia_3_300()

saveRDS(my_sim_busia_out_3_300, file = "analysis/IVER_ONE/output/out_busia_3_300.rds")

#3x600
ento_mod_busia_3_600 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
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

my_sim_busia_3_600 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_busia_list, ento_mod_busia_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_busia)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "3_600"))}, simplify = F))
  return(out_df_antag)

}

my_sim_busia_out_3_600 <- my_sim_busia_3_600()

saveRDS(my_sim_busia_out_3_600, file = "analysis/IVER_ONE/output/out_busia_3_600.rds")

ggplot(my_sim_busia_out_3_300, aes(x = t, y = EIRout))+
  geom_line()+
  geom_line(data = my_sim_busia_out_3_600, aes(x = t, y = EIRout), col = "red")+
  geom_line(data = my_sim_busia_out_1_400, aes(x = t , y = EIRout), col = "blue")+
  geom_line(data = my_sim_busia_out_1_800, aes(x = t , y = EIRout), col = "green")

#HOMA BAY####

#1x400
ento_mod_hb_1_400 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
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

my_sim_hb_1_400 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_homa_bay_list, ento_mod_hb_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_homa_bay)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA_BAY", int = "1_400"))}, simplify = F))
  return(out_df_antag)

}

my_sim_hb_out_1_400 <- my_sim_hb_1_400()

saveRDS(my_sim_hb_out_1_400, file = "analysis/IVER_ONE/output/out_hb_1_400.rds")

#1x800
ento_mod_hb_1_800 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
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

my_sim_hb_1_800 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_homa_bay_list, ento_mod_hb_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_homa_bay)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA_BAY", int = "1_800"))}, simplify = F))
  return(out_df_antag)

}

my_sim_hb_out_1_800 <- my_sim_hb_1_800()

saveRDS(my_sim_hb_out_1_800, file = "analysis/IVER_ONE/output/out_hb_1_800.rds")

#3x300
ento_mod_hb_3_300 <- function(itn_type_ivm_param){
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
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

my_sim_hb_3_300 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_homa_bay_list, ento_mod_hb_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_homa_bay)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA_BAY", int = "3_300"))}, simplify = F))
  return(out_df_antag)

}

my_sim_hb_out_3_300 <- my_sim_hb_3_300()

saveRDS(my_sim_hb_out_3_300, file = "analysis/IVER_ONE/output/out_hb_3_300.rds")

#3x600
ento_mod_hb_3_600 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
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

my_sim_hb_3_600 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_homa_bay_list, ento_mod_hb_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_homa_bay)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA_BAY", int = "3_600"))}, simplify = F))
  return(out_df_antag)

}

my_sim_hb_out_3_600 <- my_sim_hb_3_600()

saveRDS(my_sim_hb_out_3_600, file = "analysis/IVER_ONE/output/out_hb_3_600.rds")

ggplot(my_sim_hb_out_3_600, aes(x = t, y = EIRout))+
  geom_line()+
  geom_line(data = my_sim_hb_out_3_300, aes(x = t, y = EIRout), col = "red")+
  geom_line(data = my_sim_hb_out_1_400, aes(x = t, y = EIRout), col = "green")+
  geom_line(data = my_sim_hb_out_1_800, aes(x = t, y = EIRout), col = "blue")


#MIGORI####

#1x400
ento_mod_migori_1_400 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
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

my_sim_migori_1_400 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_migori_list, ento_mod_migori_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_migori)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "1_400"))}, simplify = F))
  return(out_df_antag)

}

my_sim_migori_out_1_400 <- my_sim_migori_1_400()

saveRDS(my_sim_migori_out_1_400, file = "analysis/IVER_ONE/output/out_migori_1_400.rds")

#1x800
ento_mod_migori_1_800 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
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

my_sim_migori_1_800 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_migori_list, ento_mod_migori_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_migori)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "1_800"))}, simplify = F))
  return(out_df_antag)

}

my_sim_migori_out_1_800 <- my_sim_migori_1_800()

saveRDS(my_sim_migori_out_1_800, file = "analysis/IVER_ONE/output/out_migori_1_800.rds")

#3x300
ento_mod_migori_3_300 <- function(itn_type_ivm_param){
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
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

my_sim_migori_3_300 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_migori_list, ento_mod_migori_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_migori)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "3_300"))}, simplify = F))
  return(out_df_antag)

}

my_sim_migori_out_3_300 <- my_sim_migori_3_300()

saveRDS(my_sim_migori_out_3_300, file = "analysis/IVER_ONE/output/out_migori_3_300.rds")

#3x600
ento_mod_migori_3_600 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
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

my_sim_migori_3_600 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_migori_list, ento_mod_migori_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_migori)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "3_600"))}, simplify = F))
  return(out_df_antag)

}

my_sim_migori_out_3_600 <- my_sim_migori_3_600()

saveRDS(my_sim_migori_out_3_600, file = "analysis/IVER_ONE/output/out_migori_3_600.rds")


#SIAYA####
ento_mod_siaya_1_400 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
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

my_sim_siaya_1_400 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_siaya_list, ento_mod_siaya_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_siaya)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "1_400"))}, simplify = F))
  return(out_df_antag)

}

my_sim_siaya_out_1_400 <- my_sim_siaya_1_400()

saveRDS(my_sim_siaya_out_1_400, file = "analysis/IVER_ONE/output/out_siaya_1_400.rds")

#1x800
ento_mod_siaya_1_800 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
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

my_sim_siaya_1_800 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_siaya_list, ento_mod_siaya_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_siaya)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "1_800"))}, simplify = F))
  return(out_df_antag)

}

my_sim_siaya_out_1_800 <- my_sim_siaya_1_800()

saveRDS(my_sim_siaya_out_1_800, file = "analysis/IVER_ONE/output/out_siaya_1_800.rds")

#3x300
ento_mod_siaya_3_300 <- function(itn_type_ivm_param){
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
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

my_sim_siaya_3_300 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_siaya_list, ento_mod_siaya_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_siaya)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "3_300"))}, simplify = F))
  return(out_df_antag)

}

my_sim_siaya_out_3_300 <- my_sim_siaya_3_300()

saveRDS(my_sim_siaya_out_3_300, file = "analysis/IVER_ONE/output/out_siaya_3_300.rds")

#3x600
ento_mod_siaya_3_600 <- function(itn_type_ivm_param){
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
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
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

my_sim_siaya_3_600 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(mod_siaya_list, ento_mod_siaya_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_df_siaya)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "3_600"))}, simplify = F))
  return(out_df_antag)

}
Sys.time()
my_sim_siaya_out_3_600 <- my_sim_siaya_3_600()
Sys.time()

saveRDS(my_sim_siaya_out_3_600, file = "analysis/IVER_ONE/output/out_siaya_3_600.rds")

ggplot(my_sim_migori_out_3_600, aes(x = t, y = EIRout))+
  geom_line()+
  geom_line(data = my_sim_migori_out_3_300, aes(x = t, y = EIRout), col = "red")+
  geom_line(data = my_sim_migori_out_1_400, aes(x = t, y = EIRout), col = "green")+
  geom_line(data = my_sim_migori_out_1_800, aes(x = t, y = EIRout), col = "blue")

ggplot(my_sim_siaya_out_3_600, aes(x = t, y = EIRout))+
  geom_line()+
  geom_line(data = my_sim_siaya_out_3_300, aes(x = t, y = EIRout), col = "red")+
  geom_line(data = my_sim_siaya_out_1_400, aes(x = t, y = EIRout), col = "green")+
  geom_line(data = my_sim_siaya_out_1_800, aes(x = t, y = EIRout), col = "blue")


ggplot(my_sim_busia_out_3_600, aes(x = t, y = EIRout))+
  geom_line()+
  geom_line(data = my_sim_busia_out_3_300, aes(x = t, y = EIRout), col = "red")+
  geom_line(data = my_sim_busia_out_1_400, aes(x = t, y = EIRout), col = "green")+
  geom_line(data = my_sim_busia_out_1_800, aes(x = t, y = EIRout), col = "blue")


#run exp decay model
wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")


ivm_params_kwale_df <- dplyr::cross_join(mod_df_kwale, malsim_odin)
ivm_params_busia_df <- dplyr::cross_join(mod_df_busia, malsim_odin)
ivm_params_hb_df <- dplyr::cross_join(mod_df_homa_bay, malsim_odin)
ivm_params_migori_df <-dplyr::cross_join(mod_df_migori, malsim_odin)
ivm_params_siaya_df <-dplyr::cross_join(mod_df_siaya, malsim_odin)


#make lists
ivm_params_kwale_list <- list()
for(i in seq_len(nrow(ivm_params_kwale_df))){
  ivm_params_kwale_list[[i]] <- as.numeric(ivm_params_kwale_df[i,])
}


ivm_params_busia_list <- list()
for(i in seq_len(nrow(ivm_params_busia_df))){
  ivm_params_busia_list[[i]] <- as.numeric(ivm_params_busia_df[i,])
}

ivm_params_hb_list <- list()
for(i in seq_len(nrow(ivm_params_hb_df))){
  ivm_params_hb_list[[i]] <- as.numeric(ivm_params_hb_df[i,])
}

ivm_params_migori_list <- list()
for(i in seq_len(nrow(ivm_params_migori_df))){
  ivm_params_migori_list[[i]] <- as.numeric(ivm_params_migori_df[i,])
}

ivm_params_siaya_list <- list()
for(i in seq_len(nrow(ivm_params_siaya_df))){
  ivm_params_siaya_list[[i]] <- as.numeric(ivm_params_siaya_df[i,])
}

#KWALE EXP DECAY RUNS

#1x400
exp_decay_ento_mod_kwale_1_400 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Kwale", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_kwale_1_400_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_kwale_list, exp_decay_ento_mod_kwale_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_kwale_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "1_400_exp_decay"))}, simplify = F))
  return(out_df_antag)

}


#1x800
exp_decay_ento_mod_kwale_1_800 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Kwale", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_kwale_1_800_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_kwale_list, exp_decay_ento_mod_kwale_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_kwale_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "1_800_exp_decay"))}, simplify = F))
  return(out_df_antag)

}
#3x300
exp_decay_ento_mod_kwale_3_300 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_kwale_3_300_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_kwale_list, exp_decay_ento_mod_kwale_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_kwale_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "3_300_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x600
exp_decay_ento_mod_kwale_3_600 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Kwale", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_kwale_3_600_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_kwale_list, exp_decay_ento_mod_kwale_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_kwale_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "KWALE", int = "3_600_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#BUSIA

#1x400
exp_decay_ento_mod_busia_1_400 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Busia", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_busia_1_400_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_busia_list, exp_decay_ento_mod_busia_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_busia_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "1_400_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#1x800
exp_decay_ento_mod_busia_1_800 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Busia", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_busia_1_800_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_busia_list, exp_decay_ento_mod_busia_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_busia_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "1_800_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x300
exp_decay_ento_mod_busia_3_300 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_busia_3_300_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_busia_list, exp_decay_ento_mod_busia_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_busia_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "3_300_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x600
exp_decay_ento_mod_busia_3_600 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Busia", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_busia_3_600_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_busia_list, exp_decay_ento_mod_busia_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_busia_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "BUSIA", int = "3_600_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#HOMA BAY####
#1x400

exp_decay_ento_mod_hb_1_400 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Homa Bay", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_hb_1_400_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_hb_list, exp_decay_ento_mod_hb_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_hb_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA BAY", int = "1_400_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#1x800
exp_decay_ento_mod_hb_1_800 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Homa Bay", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_hb_1_800_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_hb_list, exp_decay_ento_mod_hb_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_hb_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA BAY", int = "1_800_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x300
exp_decay_ento_mod_hb_3_300 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_hb_3_300_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_hb_list, exp_decay_ento_mod_hb_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_hb_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA BAY", int = "3_300_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x600
exp_decay_ento_mod_hb_3_600 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Homa Bay", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_hb_3_600_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_hb_list, exp_decay_ento_mod_hb_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_hb_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "HOMA BAY", int = "3_600_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#MIGORI

#1x400
exp_decay_ento_mod_migori_1_400 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Migori", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_migori_1_400_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_migori_list, exp_decay_ento_mod_migori_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_migori_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "1_400_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#1x800
exp_decay_ento_mod_migori_1_800 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Migori", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_migori_1_800_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_migori_list, exp_decay_ento_mod_migori_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_migori_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "1_800_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x300
exp_decay_ento_mod_migori_3_300 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_migori_3_300_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_migori_list, exp_decay_ento_mod_migori_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_migori_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "3_300_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x600
exp_decay_ento_mod_migori_3_600 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Migori", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_migori_3_600_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_migori_list, exp_decay_ento_mod_migori_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_migori_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "MIGORI", int = "3_600_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#SIAYA####
#1x400
exp_decay_ento_mod_siaya_1_400 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Siaya", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_400$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_400$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_400$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_400$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_siaya_1_400_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_siaya_list, exp_decay_ento_mod_siaya_1_400)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_siaya_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "1_400_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#1x800
exp_decay_ento_mod_siaya_1_800 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Siaya", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_800$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_800$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_800$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_800$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_siaya_1_800_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_siaya_list, exp_decay_ento_mod_siaya_1_800)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_siaya_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "1_800_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x300
exp_decay_ento_mod_siaya_3_300 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
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
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_300$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_siaya_3_300_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_siaya_list, exp_decay_ento_mod_siaya_3_300)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_siaya_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "3_300_exp_decay"))}, simplify = F))
  return(out_df_antag)

}

#3x600
exp_decay_ento_mod_siaya_3_600 <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    country = "Kenya", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Siaya", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1_3_600$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1_3_600$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1_3_600$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    IVRM_start = ivm_parms1_3_600$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in,
    ivm_min_age = 5, #redundant param
    ivm_max_age = 90 #redundant param
  )
  return(output)
}

my_sim_siaya_3_600_exp_decay <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_params_siaya_list, exp_decay_ento_mod_siaya_3_600)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_siaya_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIRout, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, ivm_cov, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "SIAYA", int = "3_600_exp_decay"))}, simplify = F))
  return(out_df_antag)

}






#RUN MODELS

#kwale

Sys.time()

my_sim_kwale_out_1_400_exp_decay <- my_sim_kwale_1_400_exp_decay()

Sys.time()


saveRDS(my_sim_kwale_out_1_400_exp_decay, file = "analysis/IVER_ONE/output/out_kwale_1_400_exp_decay.rds")


my_sim_kwale_out_1_800_exp_decay <- my_sim_kwale_1_800_exp_decay()
saveRDS(my_sim_kwale_out_1_800_exp_decay, file = "analysis/IVER_ONE/output/out_kwale_1_800_exp_decay.rds")

my_sim_kwale_out_3_300_exp_decay <- my_sim_kwale_3_300_exp_decay()
saveRDS(my_sim_kwale_out_3_300_exp_decay, file = "analysis/IVER_ONE/output/out_kwale_3_300_exp_decay.rds")

my_sim_kwale_out_3_600_exp_decay <- my_sim_kwale_3_600_exp_decay()
saveRDS(my_sim_kwale_out_3_600_exp_decay, file = "analysis/IVER_ONE/output/out_kwale_3_600_exp_decay.rds")

#BUSIA
my_sim_busia_out_1_400_exp_decay <- my_sim_busia_1_400_exp_decay()
saveRDS(my_sim_busia_out_1_400_exp_decay, file = "analysis/IVER_ONE/output/out_busia_1_400_exp_decay.rds")

my_sim_busia_out_1_800_exp_decay <- my_sim_busia_1_800_exp_decay()
saveRDS(my_sim_busia_out_1_800_exp_decay, file = "analysis/IVER_ONE/output/out_busia_1_800_exp_decay.rds")

my_sim_busia_out_3_300_exp_decay <- my_sim_busia_3_300_exp_decay()
saveRDS(my_sim_busia_out_3_300_exp_decay, file = "analysis/IVER_ONE/output/out_busia_3_300_exp_decay.rds")

my_sim_busia_out_3_600_exp_decay <- my_sim_busia_3_600_exp_decay()
saveRDS(my_sim_busia_out_3_600_exp_decay, file = "analysis/IVER_ONE/output/out_busia_3_600_exp_decay.rds")

#HOMA BAY
my_sim_hb_out_1_400_exp_decay <- my_sim_hb_1_400_exp_decay()
saveRDS(my_sim_hb_out_1_400_exp_decay, file = "analysis/IVER_ONE/output/out_hb_1_400_exp_decay.rds")

my_sim_hb_out_1_800_exp_decay <- my_sim_hb_1_800_exp_decay()
saveRDS(my_sim_hb_out_1_800_exp_decay, file = "analysis/IVER_ONE/output/out_hb_1_800_exp_decay.rds")

my_sim_hb_out_3_300_exp_decay <- my_sim_hb_3_300_exp_decay()
saveRDS(my_sim_hb_out_3_300_exp_decay, file = "analysis/IVER_ONE/output/out_hb_3_300_exp_decay.rds")

my_sim_hb_out_3_600_exp_decay <- my_sim_hb_3_600_exp_decay()
saveRDS(my_sim_hb_out_3_600_exp_decay, file = "analysis/IVER_ONE/output/out_hb_3_600_exp_decay.rds")

#MIGORI
#1x400
my_sim_migori_out_1_400_exp_decay <- my_sim_migori_1_400_exp_decay()
saveRDS(my_sim_migori_out_1_400_exp_decay, file = "analysis/IVER_ONE/output/out_migori_1_400_exp_decay.rds")

#1x800
my_sim_migori_out_1_800_exp_decay <- my_sim_migori_1_800_exp_decay()
saveRDS(my_sim_migori_out_1_800_exp_decay, file = "analysis/IVER_ONE/output/out_migori_1_800_exp_decay.rds")

#3x300
my_sim_migori_out_3_300_exp_decay <- my_sim_migori_3_300_exp_decay()
saveRDS(my_sim_migori_out_3_300_exp_decay, file = "analysis/IVER_ONE/output/out_migori_3_300_exp_decay.rds")

#3x600
my_sim_migori_out_3_600_exp_decay <- my_sim_migori_3_600_exp_decay()
saveRDS(my_sim_migori_out_3_600_exp_decay, file = "analysis/IVER_ONE/output/out_migori_3_600_exp_decay.rds")

#SIAYA###
#1x400
my_sim_siaya_out_1_400_exp_decay <- my_sim_siaya_1_400_exp_decay()
saveRDS(my_sim_siaya_out_1_400_exp_decay, file = "analysis/IVER_ONE/output/out_siaya_1_400_exp_decay.rds")

#1x800
my_sim_siaya_out_1_800_exp_decay <- my_sim_siaya_1_800_exp_decay()
saveRDS(my_sim_siaya_out_1_800_exp_decay, file = "analysis/IVER_ONE/output/out_siaya_1_800_exp_decay.rds")

#3x300
my_sim_siaya_out_3_300_exp_decay <- my_sim_siaya_3_300_exp_decay()
saveRDS(my_sim_siaya_out_3_300_exp_decay, file = "analysis/IVER_ONE/output/out_siaya_3_300_exp_decay.rds")

#3x600
my_sim_siaya_out_3_600_exp_decay <- my_sim_siaya_3_600_exp_decay()
saveRDS(my_sim_siaya_out_3_600_exp_decay, file = "analysis/IVER_ONE/output/out_siaya_3_600_exp_decay.rds")

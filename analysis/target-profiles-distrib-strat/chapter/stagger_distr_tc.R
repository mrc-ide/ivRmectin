#model staggered distribution by altering hazard ratios.

#we assume equal coverage at each distribution

#modify the hazard ratio curve
require(tidyverse)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")

stag_10 <- 11
times_10_d <- seq(1, stag_10, length.out = 3)

ivm_haz <- ivm_haz %>%
  select(-IVM_400_1_HS) %>%
  filter(between(Day, 1, 23))


shift_days_10_MDA_3 <- times_10_d[3] - times_10_d[1]

max_day <- max(ivm_haz$Day)
new_days <- (max_day + 1):(max_day + shift_days_10_MDA_3)

extra_rows_10d <- data.frame(Day = new_days, IVM_300_3_HS = 1)

df_extended_10d <- rbind(ivm_haz, extra_rows_10d)

#create new shifted column - group 2 - and fill with HR 1 first, impact impact

df_extended_10d$group3 <- 1 # HR1 - no impact

df_extended_10d$group3[(shift_days_10_MDA_3 + 1):nrow(df_extended_10d)] <- df_extended_10d$IVM_300_3_HS[1:(nrow(df_extended_10d) - shift_days_10_MDA_3)]

shift_days_10_MDA_2 <- times_10_d[2] - times_10_d[1]

df_extended_10d$group2 <- 1
df_extended_10d$group2[(shift_days_10_MDA_2 + 1):nrow(df_extended_10d)] <- df_extended_10d$IVM_300_3_HS[1:(nrow(df_extended_10d) - shift_days_10_MDA_2)]

df_extended_10d <- df_extended_10d %>%
  mutate(group1 = IVM_300_3_HS) %>%
  select(Day, IVM_300_3_HS, group1, group2, group3)

df_extended_10d_new_HR <- df_extended_10d %>%
  #mutate(across(everything(), ~ replace_na(.x, 1))) %>% #no killing effect of others as MDA not started yet
  rowwise() %>%
  mutate(HR_use = mean(c(group1, group2, group3))) %>%
  mutate(stagger = "10d")

#for the 30d strategy

stag_30 <- 31
times_30_d <- seq(1, stag_30, length.out = 3)

shift_days_30_MDA_3 <- times_30_d[3] - times_30_d[1]

max_day <- max(ivm_haz$Day)
new_days <- (max_day + 1):(max_day + shift_days_30_MDA_3)

extra_rows_30d <- data.frame(Day = new_days, IVM_300_3_HS = 1)

df_extended_30d <- rbind(ivm_haz, extra_rows_30d)

#create new shifted column - group 2 - and fill with HR 1 first, no impact

df_extended_30d$group3 <- 1 #HR 1 - no impact

df_extended_30d$group3[(shift_days_30_MDA_3 + 1):nrow(df_extended_30d)] <- df_extended_30d$IVM_300_3_HS[1:(nrow(df_extended_30d) - shift_days_30_MDA_3)]

shift_days_30_MDA_2 <- times_30_d[2] - times_30_d[1]

df_extended_30d$group2 <- 1
df_extended_30d$group2[(shift_days_30_MDA_2 + 1):nrow(df_extended_30d)] <- df_extended_30d$IVM_300_3_HS[1:(nrow(df_extended_30d) - shift_days_30_MDA_2)]

df_extended_30d <- df_extended_30d %>%
  mutate(group1 = IVM_300_3_HS) %>%
  select(Day, IVM_300_3_HS, group1, group2, group3)

df_extended_30d_new_HR <- df_extended_30d %>%
  #mutate(across(everything(), ~ replace_na(.x, 1))) %>% #so can take average
  rowwise() %>%
  mutate(HR_use = mean(c(group1, group2, group3))) %>%
  mutate(stagger = "30d")

HR_all <- rbind(df_extended_10d_new_HR, df_extended_30d_new_HR)

ggplot(HR_all, aes(x = Day, y = HR_use, col = as.factor(stagger)))+
  geom_line()+
  theme_bw()

#we now pass this into the model

time_period <- 365*15 #long run to get to eqm
mda_int <- 30
ivm_cov_in = c(0.1, 0.9)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

#itn_on <- 100 #introduce nets 100 days into simulation

#modelling with extreme density dependence

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

#distr plan
start <- (365*5)+200 #time for starting distr, help get to reach eqm, and to distr at right time when seasonality introduced

Q0_in <- c(0.21, 0.71, 0.92, 0.94)
init_EIR_in <- c(2, 100)
target_cov <- 0.9

df_var_all <- expand.grid(Q0 = Q0_in, ivm_cov = target_cov, init_EIR = init_EIR_in)
#df <- df[1,]
my_list_all <- list()
for (i in seq_len(nrow(df_var_all))){
  my_list_all[[i]] <- as.numeric(df_var_all[i,])
}



ivm_parms_10d_stag <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  time_period = time_period,
  hazard_profile = df_extended_10d_new_HR$HR_use,
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

mod_10d_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_10d_stag$ttt,
    eff_len = ivm_parms_10d_stag$eff_len,
    haz = ivm_parms_10d_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_10d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_10d_stag$ivm_max_age,
    IVRM_start = ivm_parms_10d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_10d_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_10d_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "10d-stagger"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_10d_stag <- my_sim_mod_10d_stag()

#then 30d stagger
ivm_parms_30d_stag <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  time_period = time_period,
  hazard_profile = df_extended_30d_new_HR$HR_use,
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

mod_30d_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_30d_stag$ttt,
    eff_len = ivm_parms_30d_stag$eff_len,
    haz = ivm_parms_30d_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_30d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_30d_stag$ivm_max_age,
    IVRM_start = ivm_parms_30d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_30d_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_30d_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "30d-stagger"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_30d_stag <- my_sim_mod_30d_stag()

#all-in-one
ivm_parms_all <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

mod_all <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_all <- function(){
  mod_out_list <- lapply(my_list_all, mod_all)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "all-in-one"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_all <- my_sim_mod_all()

#baseline scenario

mod_baseline <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = 0, #make coverage 0
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_baseline <- function(){
  mod_out_list <- lapply(my_list_all, mod_baseline)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "baseline"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_baseline <- my_sim_mod_baseline()

df_mod_distr <- do.call("rbind", list(df_mod_10d_stag, df_mod_30d_stag, df_mod_all, df_mod_baseline))

write_rds(df_mod_distr, file = "analysis/target-profiles-distrib-strat/chapter/output/df_distr_HR.rds")

#with seasonality

#10d Stagger

mod_10d_stag_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_10d_stag$ttt,
    eff_len = ivm_parms_10d_stag$eff_len,
    haz = ivm_parms_10d_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_10d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_10d_stag$ivm_max_age,
    IVRM_start = ivm_parms_10d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_10d_stag_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_10d_stag_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "10d-stagger-Sen"))}, simplify = F))
  return(mod_df)
}

df_mod_10d_stag_Sen <- my_sim_mod_10d_stag_Sen()

#over 30 days

mod_30d_stag_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_30d_stag$ttt,
    eff_len = ivm_parms_30d_stag$eff_len,
    haz = ivm_parms_30d_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_30d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_30d_stag$ivm_max_age,
    IVRM_start = ivm_parms_30d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_30d_stag_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_30d_stag_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "30d-stagger-Sen"))}, simplify = F))
  return(mod_df)
}

df_mod_30d_stag_Sen <- my_sim_mod_30d_stag_Sen()

#all in one, with seasonality
mod_all_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_all_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_all_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "all-in-one-Sen"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_all_Sen <- my_sim_mod_all_Sen()

#baseline
mod_baseline_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = 0, #make coverage 0
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_baseline_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_baseline_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "baseline-Sen"))}, simplify = F))
  return(mod_df)
}

df_mod_baseline_Sen <- my_sim_mod_baseline_Sen()

df_mod_distr_Sen <- do.call("rbind", list(df_mod_10d_stag_Sen, df_mod_30d_stag_Sen, df_mod_all_Sen, df_mod_baseline_Sen))

write_rds(df_mod_distr_Sen, file = "analysis/target-profiles-distrib-strat/chapter/output/df_distr_HR_Sen.rds")

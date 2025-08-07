#looking at how best to kill mosquitoes

#to kill 100 mosquitoes, should we do it in 1 day or over a longer period of time?

#first endectocide has a killing effect for 23 days, and increases the death rate from 0.1 to 0.2 (for example)
#so endec_mu for endectocide 1 is 0.1

#other products have shorter or longer periods of killing, what endec_mu is required to match the number of dead mosquitoes from product 1
#with the following killing times:
#1 7 days
#2 14 days
#3 30 days
#4 6 months (180 days)

#script to compare different model types
devtools::load_all()
require(tidyverse)
#how much does mu_h change across different values of Q0 and ivm_cov?

# Provide a value of the annual EIR for this model run
#init_EIR_vec <- c(2, 25, 100) #low - 2, moderate - 15, high - 120 --> Ellie: low = 2, med = 25, high = 100

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years, turn ivermectin on when nets are 6m, 1y, 2.5yo
time_period <- 365*8
mda_int <- 30
ivm_cov_in = c(0.1, 0.9)

#ivm on when nets are 6 months old
net_seq <- seq(100, 3650, by = 3*365)
IVM_begin1 <- net_seq[3]+180 # 6 months into new net distribution
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

#itn_on <- 100 #introduce nets 100 days into simulation


runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#here the hazard profile is 1 everyday

#set up the different start times (coverage is redundant here, we are just operating a switch)

ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, 23),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms2 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, 7),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms3 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, 14),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms4 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, 30),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms5 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, 180),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

#for the first model (23 day killing), let's take an endec_mu of 0.01 --> endectocide mortality rate is 0.01 + 0.132 = 0.142
#and mosquitoes live on for average, 7.04 days

mod1 <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 25,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms1$ttt,
    eff_len = ivm_parms1$eff_len,
    haz = ivm_parms1$haz,
    ivm_cov_par = ivm_parms1$ivm_cov_par,
    ivm_min_age = ivm_parms1$ivm_min_age,
    ivm_max_age = ivm_parms1$ivm_max_age,
    IVRM_start = ivm_parms1$IVRM_start,
    Q0 = 0.9,
    endec_mu = 0.1)

mod1_df <- as.data.frame(runfun(mod1)) %>%
  select(t, mu, mv, mv_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr)

ggplot(mod1_df, aes(x = t, y = mv))+
  geom_line() #got the drops

mod1_endec_killed <- mod1_df %>%
  filter(between(t, IVM_start[1], IVM_start[1]+23)) %>%
  select(t, mv_dead)

ggplot(mod1_endec_killed, aes(x = t, y = mv_dead))+
  geom_line()

#then for endectocides that kill for 7 days, pass in a range of endec_mu
endec_mu_in <- seq(0, 1, 0.01)
df_var <- data.frame(endec_mu_in = endec_mu_in)
my_list <- list()

for (i in seq_len(nrow(df_var))){
  my_list[[i]] <- as.numeric(df_var[i,])
}

mod2 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 25,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms2$ttt,
    eff_len = ivm_parms2$eff_len,
    haz = ivm_parms2$haz,
    ivm_cov_par = ivm_parms2$ivm_cov_par,
    ivm_min_age = ivm_parms2$ivm_min_age,
    ivm_max_age = ivm_parms2$ivm_max_age,
    IVRM_start = ivm_parms2$IVRM_start,
    Q0 = 0.9,
    endec_mu = endec_mu_in
  )
  return(output)
}

my_sim_mod2 <- function(){
  mod2_out_list <- lapply(my_list, mod2)
  res_mod2_out <- lapply(mod2_out_list, runfun)
  mod2_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod2_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mv_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu", eff_len = "7"))}, simplify = F))
  return(mod2_df)
}

df_mod2 <- my_sim_mod2()

#df_mod2_dead <- df_mod2 %>%
#  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
#  select(t, mv_dead, ref)
#
#mod2_list <- split(df_mod2_dead, f = df_mod2_dead$ref)
#error <- numeric()
#for (i in 1:length(endec_mu_in)){
#  error <- c(error, sum(mod1_endec_killed$mv_dead - mod2_list[[i]]$mv_dead)^2)
#}
#index <- which.min(error)
#endec_mu_in[index] #so an endec_mu of 0.62 over 7 days is required to match endec_mu of 0.01 over 23 days
#range(df_mod2_dead$ref)


#then repeat this slog for all the others

#14 days
mod3 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 25,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms3$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = ivm_parms3$ivm_cov_par,
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start,
    Q0 = 0.9,
    endec_mu = endec_mu_in
  )
  return(output)
}

my_sim_mod3 <- function(){
  mod3_out_list <- lapply(my_list, mod3)
  res_mod3_out <- lapply(mod3_out_list, runfun)
  mod3_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod3_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mv_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu", eff_len = "14"))}, simplify = F))
  return(mod3_df)
}

df_mod3 <- my_sim_mod3()

#fitting
#df_mod3_dead <- df_mod3 %>%
#  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
#  select(t, mv_dead, ref)
#
#mod3_list <- split(df_mod3_dead, f = df_mod3_dead$ref)
#error <- numeric()
#for (i in 1:length(endec_mu_in)){
#  error <- c(error, sum(mod1_endec_killed$mv_dead - mod3_list[[i]]$mv_dead)^2)
#}
#index <- which.min(error)
#endec_mu_in[index] #so an endec_mu of 0.19 over 14 days is required to match endec_mu of 0.01 over 23 days


#30 days
mod4 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 25,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms4$ttt,
    eff_len = ivm_parms4$eff_len,
    haz = ivm_parms4$haz,
    ivm_cov_par = ivm_parms4$ivm_cov_par,
    ivm_min_age = ivm_parms4$ivm_min_age,
    ivm_max_age = ivm_parms4$ivm_max_age,
    IVRM_start = ivm_parms4$IVRM_start,
    Q0 = 0.9,
    endec_mu = endec_mu_in
  )
  return(output)
}

my_sim_mod4 <- function(){
  mod4_out_list <- lapply(my_list, mod4)
  res_mod4_out <- lapply(mod4_out_list, runfun)
  mod4_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod4_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mv_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu", eff_len = "30"))}, simplify = F))
  return(mod4_df)
}

df_mod4 <- my_sim_mod4()

#fitting
#df_mod4_dead <- df_mod4 %>%
#  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
#  select(t, mv_dead, ref)
#
#mod4_list <- split(df_mod4_dead, f = df_mod4_dead$ref)
#error <- numeric()
#for (i in 1:length(endec_mu_in)){
#  error <- c(error, sum(mod1_endec_killed$mv_dead - mod4_list[[i]]$mv_dead)^2)
#}
#index <- which.min(error)
#endec_mu_in[index] #so an endec_mu of 0.08 over 30 days is required to match endec_mu of 0.01 over 23 days


#180 days
mod5 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 25,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms5$ttt,
    eff_len = ivm_parms5$eff_len,
    haz = ivm_parms5$haz,
    ivm_cov_par = ivm_parms5$ivm_cov_par,
    ivm_min_age = ivm_parms5$ivm_min_age,
    ivm_max_age = ivm_parms5$ivm_max_age,
    IVRM_start = ivm_parms5$IVRM_start,
    Q0 = 0.9,
    endec_mu = endec_mu_in
  )
  return(output)
}

my_sim_mod5 <- function(){
  mod5_out_list <- lapply(my_list, mod5)
  res_mod5_out <- lapply(mod5_out_list, runfun)
  mod5_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod5_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mv_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu", eff_len = "180"))}, simplify = F))
  return(mod5_df)
}

df_mod5 <- my_sim_mod5()

#fitting
#df_mod5_dead <- df_mod5 %>%
#  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
#  select(t, mv_dead, ref)
#
#mod5_list <- split(df_mod5_dead, f = df_mod5_dead$ref)
#error <- numeric()
#for (i in 1:length(endec_mu_in)){
#  error <- c(error, sum(mod1_endec_killed$mv_dead - mod5_list[[i]]$mv_dead)^2)
#}
#index <- which.min(error)
#endec_mu_in[index] #so an endec_mu of 0.08 over 30 days is required to match endec_mu of 0.01 over 23 days


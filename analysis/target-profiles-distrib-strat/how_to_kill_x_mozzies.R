#analysis to work out whether it is best to kill the same number of mosquitoes in a short period of time or long period of time


#looking at how best to kill mosquitoes

#to kill 100 mosquitoes, should we do it in 1 day or over a longer period of time to have biggest impact on transmission (EIR)?

#first endectocide has a killing effect for 23 days, and increases the death rate from 0.1 to 0.2 (for example)
#so endec_mu for endectocide 1 is 0.1

#other products have shorter or longer periods of killing, what endec_mu is required to match the number of dead mosquitoes from product 1
#with the following killing times:
#1 23 days
#2 20 days
#3 26 days
#4 40 days

#script to compare different model types
devtools::load_all()
require(tidyverse)
#how much does mu_h change across different values of Q0 and ivm_cov?

# Provide a value of the annual EIR for this model run
#init_EIR_vec <- c(2, 25, 100) #low - 2, moderate - 15, high - 120 --> Ellie: low = 2, med = 25, high = 100

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years, turn ivermectin on when nets are 6m, 1y, 2.5yo
time_period <- 365*3
mda_int <- 30
#ivm_cov_in = c(0.1, 0.9)

#ivm on when nets are 6 months old

IVM_begin1 <- 500
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

#itn_on <- 100 #introduce nets 100 days into simulation


runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

#here the hazard profile is 1 everyday
#just doing MDA for one month (not every month for 3 months)

#set up the different start times (coverage is redundant here, we are just operating a switch)
eff_len_endec1 <- 23
ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, eff_len_endec1),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)


eff_len_endec2 <- 20
ivm_parms2 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, eff_len_endec2),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)


eff_len_endec3 <- 26
ivm_parms3 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, eff_len_endec3),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)


eff_len_endec4 <- 40
ivm_parms4 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, eff_len_endec4),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

#using models with density dependence on

#counterfactual: no interventions
mod0 <- ivRmectin:::create_r_model(
  odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay3.R", package = "ivRmectin"),
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
  IVRM_start = ivm_parms1$IVRM_start,
  Q0 = 0.9,
  endec_mu = 0,
  wane = 0,
  init_ft = 0)

mod0_df <- as.data.frame(runfun(mod0)) %>%
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, mvx_dead_natural, betaa, beta_larval, KL)
write_rds(mod0_df, file = "analysis/target-profiles-distrib-strat/killing_mosq_endec_mod0_exp_decay3.rds")

mod1 <- ivRmectin:::create_r_model(
  odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay4.R", package = "ivRmectin"),
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
  IVRM_start = ivm_parms1$IVRM_start,
  Q0 = 0.9,
  endec_mu = 0.09,
  wane = 0)

mod1_df <- as.data.frame(runfun(mod1)) %>%
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, Ivtot,clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural) %>%
  mutate(model_type = "endec_mu", eff_len = eff_len_endec1)


write_rds(mod1_df, file = "analysis/target-profiles-distrib-strat/killing_mosq_endec_mod1_exp_decay3.rds")



#then for endectocides that kill for 23 days, pass in a range of endec_mu
endec_mu_in <- seq(0.01, 1, 0.01)
#endec_mu_in <- 0.53
df_var <- data.frame(endec_mu_in = endec_mu_in)
my_list <- list()

for (i in seq_len(nrow(df_var))){
  my_list[[i]] <- as.numeric(df_var[i,])
}

mod2 <-  function(data_in){
  endec_mu_in <- data_in[1]
  #endec_mu_in <- 0.53
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay3.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
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
    endec_mu = endec_mu_in,
    wane = 0
    #endec_mu = 0.53
  )
  return(output)
}

my_sim_mod2 <- function(){
  mod2_out_list <- lapply(my_list, mod2)
  res_mod2_out <- lapply(mod2_out_list, runfun)
  mod2_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod2_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, Ivtot,clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu", eff_len = eff_len_endec2))}, simplify = F))
  return(mod2_df)
}

df_mod2 <- my_sim_mod2()
write_rds(df_mod2, file = "analysis/target-profiles-distrib-strat/killing_mosq_endec_mod2_exp_decay3.rds")

endec_mu2 <- 0.09*23/eff_len_endec2
#endec_mu_in <- 0.53
df_var <- data.frame(endec_mu_in = endec_mu2)
my_list <- list()

for (i in seq_len(nrow(df_var))){
  my_list[[i]] <- as.numeric(df_var[i,])
}

mod2_AG <-  function(data_in){
  endec_mu_in <- data_in[1]
  #endec_mu_in <- 0.53
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay4.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
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
    endec_mu = endec_mu_in,
    wane = 0
    #endec_mu = 0.53
  )
  return(output)
}

my_sim_mod2_AG <- function(){
  mod2_out_list <- lapply(my_list, mod2_AG)
  res_mod2_out <- lapply(mod2_out_list, runfun)
  mod2_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod2_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, Ivtot,clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, model_type = "endec_mu", eff_len = eff_len_endec2))}, simplify = F))
  return(mod2_df)
}

df_mod2_AG <- my_sim_mod2_AG()



df_var <- data.frame(endec_mu_in = endec_mu_in)
my_list <- list()

for (i in seq_len(nrow(df_var))){
  my_list[[i]] <- as.numeric(df_var[i,])
}


#26 days
mod3 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay3.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
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
    wane = 0,
    endec_mu = endec_mu_in
    #endec_mu = 0.01
  )
  return(output)
}

my_sim_mod3 <- function(){
  mod3_out_list <- lapply(my_list, mod3)
  res_mod3_out <- lapply(mod3_out_list, runfun)
  mod3_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod3_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, endec_mu, EIR_tot, Ivtot,slide_prev0to5, IVRM_sr,clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, model_type = "endec_mu", eff_len = eff_len_endec3))}, simplify = F))
  return(mod3_df)
}

df_mod3 <- my_sim_mod3()
write_rds(df_mod3, file = "analysis/target-profiles-distrib-strat/killing_mosq_endec_mod3_exp_decay3.rds")


endec_mu3 <- (0.09*23)/eff_len_endec3
df_var <- data.frame(endec_mu_in = endec_mu3)
my_list <- list()

for (i in seq_len(nrow(df_var))){
  my_list[[i]] <- as.numeric(df_var[i,])
}


#26 days
mod3_AG <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay4.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
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
    wane = 0,
    endec_mu = endec_mu_in
    #endec_mu = 0.01
  )
  return(output)
}

my_sim_mod3_AG <- function(){
  mod3_out_list <- lapply(my_list, mod3_AG)
  res_mod3_out <- lapply(mod3_out_list, runfun)
  mod3_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod3_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, Ivtot,clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, model_type = "endec_mu", eff_len = eff_len_endec3))}, simplify = F))
  return(mod3_df)
}

df_mod3_AG <- my_sim_mod3_AG()


#40 days
mod4 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay3.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
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
    endec_mu = endec_mu_in,
    wane = 0
  )
  return(output)
}

my_sim_mod4 <- function(){
  mod4_out_list <- lapply(my_list, mod4)
  res_mod4_out <- lapply(mod4_out_list, runfun)
  mod4_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod4_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu", eff_len = eff_len_endec4))}, simplify = F))
  return(mod4_df)
}

df_mod4 <- my_sim_mod4()
write_rds(df_mod4, file = "analysis/target-profiles-distrib-strat/killing_mosq_endec_mod4_exp_decay3.rds")

endec_mu4 <- (0.09*23)/eff_len_endec4
df_var <- data.frame(endec_mu_in = endec_mu4)
my_list <- list()

for (i in seq_len(nrow(df_var))){
  my_list[[i]] <- as.numeric(df_var[i,])
}


#40 days
mod4_AG <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay4.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
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
    endec_mu = endec_mu_in,
    wane = 0
  )
  return(output)
}

my_sim_mod4_AG <- function(){
  mod4_out_list <- lapply(my_list, mod4_AG)
  res_mod4_out <- lapply(mod4_out_list, runfun)
  mod4_df <- do.call(rbind, sapply(1:(nrow(df_var)), function(x){
    df <- as.data.frame(res_mod4_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, Ivtot,clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, model_type = "endec_mu", eff_len = eff_len_endec4))}, simplify = F))
  return(mod4_df)
}

df_mod4_AG <- my_sim_mod4_AG()

#check fits with AG
AG_out <- do.call("rbind", list(mod1_df, df_mod2_AG, df_mod3_AG, df_mod4_AG))

AG_killed_mosq <- ggplot(AG_out, aes(x = t, y = mvx_dead, col = as.factor(eff_len)))+
  geom_line()+
  ylim(0, 80)+
  xlim(0, 700)

AG_EIR <- ggplot(AG_out, aes(x = t, y = EIR_tot, col = as.factor(eff_len)))+
  geom_line()

AG_prev <- ggplot(AG_out, aes(x = t, y = slide_prev0to5, col = as.factor(eff_len)))+
  geom_line()

AG_mv <- ggplot(AG_out, aes(x = t, y = mv, col = as.factor(eff_len)))+
  geom_line()+
  xlim(0, 700)

AG_betaa <- ggplot(AG_out, aes(x = t, y = betaa, col = as.factor(eff_len)))+
  geom_line()+
  xlim(0, 700)


cowplot::plot_grid(AG_killed_mosq, AG_EIR, AG_prev,
                   AG_mv, AG_betaa)


#compare efficacy in





#then we do fits to the number of mosquitoes killed by each endectocide.


mod1_endec_killed <- mod1_df %>%
  filter(t == IVM_start[1]+eff_len_endec1) %>%
  summarise(mvx_dead = mvx_dead, mvx_dead_natural = mvx_dead_natural)
max(mod1_endec_killed$mvx_dead) #57.4


df_mod2_dead <- df_mod2 %>%
  filter(t == IVM_start[1]+eff_len_endec2) %>%
  select(t,mv, mvx_dead, ref, endec_mu) %>%
  group_by(ref) %>%
  summarise(tot_mvx_dead = mvx_dead, endec_mu = endec_mu) #for each endec_mu, get the total number of mosquitoes killed
range(df_mod2_dead$tot_mvx_dead) #7.87 to 115
mod2_list <- split(df_mod2_dead, f = df_mod2_dead$ref)
error <- numeric()
for (i in 1:length(endec_mu_in)){
  error <- c(error, sum(mod1_endec_killed$mvx_dead - mod2_list[[i]]$tot_mvx_dead)^2)
}
index_mod2 <- which.min(error) #11
endec_mu_in[index_mod2] #0.11
range(df_mod2_dead$ref)
df_mod2_dead[index_mod2,]

best_fit_mod2 <- df_mod2 %>%
  filter(ref == index_mod2) %>%
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural) %>%
  mutate(eff_len = "20")

df_mod3_dead <- df_mod3 %>%
  filter(t == IVM_start[1]+eff_len_endec3) %>%
  select(t, mv, mvx_dead, ref, endec_mu) %>%
  group_by(ref) %>%
  summarise(mv = mv, tot_mvx_dead = mvx_dead, endec_mu)
range(df_mod3_dead$tot_mvx_dead)
mod3_list <- split(df_mod3_dead, f = df_mod3_dead$ref)
error <- numeric()
for (i in 1:length(endec_mu_in)){
  error <- c(error, sum(mod1_endec_killed$mvx_dead - mod3_list[[i]]$tot_mvx_dead)^2)
}
index_mod3 <- which.min(error) #8
endec_mu_in[index_mod3] #0.08
df_mod3_dead[index_mod3,]

df_mod3 <- readRDS("analysis/target-profiles-distrib-strat/killing_mosq_endec_mod3_exp_decay3.rds")

best_fit_mod3 <- df_mod3 %>%
  filter(ref == index_mod3) %>%
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural) %>%
  mutate(eff_len = "26")

df_mod4_dead <- df_mod4 %>%
  filter(t == IVM_start[1]+eff_len_endec4) %>%
  select(t, mv, mvx_dead, ref, endec_mu) %>%
  group_by(ref) %>%
  summarise(mv = mv, tot_mvx_dead = mvx_dead, endec_mu)
range(df_mod4_dead$tot_mvx_dead)
mod4_list <- split(df_mod4_dead, f = df_mod4_dead$ref)
error <- numeric()
for (i in 1:length(endec_mu_in)){
  error <- c(error, sum(mod1_endec_killed$mvx_dead - mod4_list[[i]]$tot_mvx_dead)^2)
}
index_mod4 <- which.min(error) #4
endec_mu_in[index_mod4] #0.04
df_mod4_dead[index_mod4,]

best_fit_mod4 <- df_mod4 %>%
  filter(ref == index_mod4) %>%
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural) %>%
  mutate(eff_len = "40")

mod1_df <- mod1_df %>%
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural) %>%
  mutate(eff_len = "23")

mod0_df <- mod0_df %>%
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, betaa, beta_larval, KL, mvx_dead_natural) %>%
  mutate(eff_len = "NA")


products <- do.call("rbind", list(mod0_df, mod1_df,
                                  best_fit_mod2, best_fit_mod3,
                                  best_fit_mod4))



mvx_dead_plot <- ggplot(products, aes(x = t, y = mvx_dead, col = as.factor(eff_len)))+
  geom_line()+
  theme_minimal()

eir_plot <- ggplot(products, aes(x = t, y = EIR_tot, col = as.factor(eff_len)))+
  geom_line()+
  theme_minimal()

mv_plot <- ggplot(products, aes(x = t, y = mv, col = as.factor(eff_len)))+
  geom_line()+
  theme_minimal()

inc_plot <- ggplot(products, aes(x = t, y = clin_inc0to5, col = as.factor(eff_len)))+
  geom_line()+
  theme_minimal()+
  ylim(0, 0.01)

prev_plot <- ggplot(products, aes(x = t, y = slide_prev0to5, col = as.factor(eff_len)))+
  geom_line()+
  theme_minimal()+
  ylim(0.5, 1)

betaa_plot <- ggplot(products, aes(x = t, y = betaa, col = as.factor(eff_len)))+
  geom_line()+
  theme_minimal()

cowplot::plot_grid(mvx_dead_plot, eir_plot, mv_plot,
                   inc_plot, prev_plot,betaa_plot, ncol = 2,
                   nrow = 3)

#TOM

#birth rate (assume to be constant)

beta<-5.7

##background mortality

mu=0.1

##population size at equibilibrium before intervention

Mv0<-beta/mu

#number of mosquitoes you need to kill per day (calculated from Andrew's equatio)

gammaT<-2

#number of mosqutioes

Mv<-seq(0,1,0.01)*Mv0

#death rate that is needed to achive gammaT

gamma<-gammaT/Mv

plot(gamma~Mv)


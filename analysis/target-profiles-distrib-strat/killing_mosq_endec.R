#looking at how best to kill mosquitoes

#to kill 100 mosquitoes, should we do it in 1 day or over a longer period of time to have biggest impact on transmission (EIR)?

#first endectocide has a killing effect for 23 days, and increases the death rate from 0.1 to 0.2 (for example)
#so endec_mu for endectocide 1 is 0.1

#other products have shorter or longer periods of killing, what endec_mu is required to match the number of dead mosquitoes from product 1
#with the following killing times:
#1 23 days
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

eff_len_endec1b <- 50
ivm_parms1b <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, eff_len_endec1b),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

eff_len_endec1c <- 50
ivm_parms1c <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start[1],
  time_period = time_period,
  hazard_profile = rep(1, eff_len_endec1b),
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



#counterfactual: no interventions
mod0 <- ivRmectin:::create_r_model(
  odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay2.R", package = "ivRmectin"),
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
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, mvx_dead_natural)
#write_rds(mod0_df, file = "analysis/target-profiles-distrib-strat/killing_mose_endec_mod0.rds")


ggplot(mod0_df, aes(x = t, y = mvx_dead_natural))+
  geom_line()

ggplot(mod0_df, aes(x = t, y = mvx_dead))+
  geom_line()

ggplot(mod0_df, aes(x = t, y = clin_inc0to5))+
  geom_line()+
  ylim(0,0.05)

#for the first model (23 day killing), let's take an endec_mu of 0.09
#and mosquitoes live on for average, 7.04 days

mod1 <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay2.R", package = "ivRmectin"),
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
  select(t, mu, mv, mvx_dead, endec_mu, EIR_tot, slide_prev0to5, IVRM_sr, clin_inc0to5, betaa, beta_larval, KL, Ivtot, eff_len,, mvx_dead_natural)


ggplot(mod1_df, aes(x = t, y = mvx_dead))+
  geom_line()

mod1_endec_killed <- mod1_df %>%
  filter(t == IVM_start[1]+eff_len_endec1) %>%
  summarise(mvx_dead = mvx_dead, mvx_dead_natural = mvx_dead_natural)


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
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay2.R", package = "ivRmectin"),
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
write_rds(df_mod2, file = "analysis/target-profiles-distrib-strat/killing_mose_endec_mod2.rds")


df_var <- data.frame(endec_mu_in = endec_mu_in)
my_list <- list()

for (i in seq_len(nrow(df_var))){
  my_list[[i]] <- as.numeric(df_var[i,])
}


#26 days
mod3 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay2.R", package = "ivRmectin"),
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
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu", eff_len = eff_len_endec3))}, simplify = F))
  return(mod3_df)
}

df_mod3 <- my_sim_mod3()
write_rds(df_mod3, file = "analysis/target-profiles-distrib-strat/killing_mose_endec_mod3.rds")

#40 days
mod4 <-  function(data_in){
  endec_mu_in <- data_in[1]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay2.R", package = "ivRmectin"),
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
write_rds(df_mod4, file = "analysis/target-profiles-distrib-strat/killing_mose_endec_mod4.rds")

#up to here

df_mod2 <- readRDS(file = "analysis/target-profiles-distrib-strat/killing_mose_endec_mod2.rds")

df_mod2_dead <- df_mod2 %>%
  filter(t == IVM_start[1]+eff_len_endec2) %>%
  select(t,mv, mvx_dead, ref, endec_mu) %>%
  group_by(ref) %>%
  summarise(tot_mvx_dead = mvx_dead, endec_mu = endec_mu) #for each endec_mu, get the total number of mosquitoes killed
range(df_mod2_dead$tot_mvx_dead) #7.87 to 130
mod2_list <- split(df_mod2_dead, f = df_mod2_dead$ref)
error <- numeric()
for (i in 1:length(endec_mu_in)){
  error <- c(error, sum(mod1_endec_killed$mvx_dead - mod2_list[[i]]$tot_mvx_dead)^2)
}
index_mod2 <- which.min(error) #11
endec_mu_in[index_mod2] #0.11
range(df_mod2_dead$ref)
df_mod2_dead[index_mod2,]

df_mod2 %>%
  filter(ref == index_mod2) %>%
  ggplot()+
  aes(x = t, y = mvx_dead)+
  geom_line()

best_fit_mod2 <- df_mod2 %>%
  filter(ref == index_mod2)

#can see that the models purge out different numbers of mosquitoes to give same number of ivermectin-killed mosquitoes
mod1_killed_mosq <- ggplot(mod1_df, aes(x = t, y = mvx_dead))+
  geom_line()+
  ggtitle("Product1: endec_mu = 0.09,eff_len = 23")+
  xlim(IVM_start[1]-10, IVM_start[1]+eff_len_endec3)


mod1_mosq <- ggplot(mod1_df, aes(x = t, y = mv))+
  geom_line()+
  ylim(0, 50)+
  ggtitle("Product1: endec_mu = 0.09,eff_len = 23")

mod1_EIR <- ggplot(mod1_df, aes(x = t, y = EIR_tot))+
  geom_line()+
  ggtitle("Product1: endec_mu = 0.09,eff_len = 23")+
  ylim(0, 200)


mod2_killed_mosq <- df_mod2 %>%
  filter(ref == index_mod2) %>%
  ggplot()+
  aes(x = t, y = mvx_dead)+
  geom_line()+
  ggtitle("Product2: endec_mu = 0.53, eff_len = 20")+
  xlim(IVM_start[1]-10, IVM_start[1]+eff_len_endec3)

mod2_mosq <- df_mod2 %>%
  filter(ref == 53) %>%
  ggplot()+
  aes(x = t, y = mv)+
  geom_line()+
  ylim(0, 50)+
  ggtitle("Product2: endec_mu = 0.53, eff_len = 20")

mod2_EIR <- df_mod2 %>%
  filter(ref == 53) %>%
  ggplot()+
  aes(x = t, y = EIR_tot)+
  geom_line()+
  ylim(0, 200)+
  ggtitle("Product2: endec_mu = 0.53, eff_len = 20")

product1_product2_plot <- cowplot::plot_grid(mod1_mosq, mod2_mosq,
                   mod1_killed_mosq, mod2_killed_mosq,
                   mod1_EIR, mod2_EIR, ncol =2,
                   nrow =3)

#impact calculation for endec 1

ggplot(mod0_df, aes(x = t, y = EIR_tot))+
  geom_line()

EIR_mod0 <- mod0_df %>%
  filter(between(t, IVM_start[1], IVM_start[3]+eff_len_endec1 )) %>%
  summarise(mean_EIR = mean(EIR_tot)) #75.58695

EIR_mod1 <- mod1_df %>%
  filter(between(t, IVM_start[1], IVM_start[3]+eff_len_endec1 )) %>%
  summarise(mean_EIR = mean(EIR_tot)) #51.23

EIR_mod2 <- df_mod2 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+eff_len_endec2 )) %>%
  summarise(mean_EIR = mean(EIR_tot)) #40.78

#EIR averted by product 1
(EIR_mod0$mean_EIR - EIR_mod1$mean_EIR)/EIR_mod0$mean_EIR #32.2% decrease
(EIR_mod0$mean_EIR - EIR_mod2$mean_EIR)/EIR_mod0$mean_EIR #46.0% decrease

#then repeat this slog for all the others


df_mod3 <- readRDS(file = "analysis/target-profiles-distrib-strat/killing_mose_endec_mod3.rds")
#fitting
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

best_fit_mod3 <- df_mod3 %>%
  filter(ref == index_mod3)

mod3_plot_mosq_killed <- df_mod3 %>%
  filter(ref == index_mod3) %>%
  ggplot()+
  aes(x = t, y = mvx_dead)+
  geom_line()+
  ggtitle("Product 3: endec_mu = 0.08, eff_len = 26")+
  xlim(IVM_start[1]-10, IVM_start[1]+eff_len_endec3)

mod3_plot_mosq <- df_mod3 %>%
  filter(ref == index_mod3) %>%
  ggplot()+
  aes(x = t, y = mv)+
  geom_line()+
  ylim(0, 50)+
  ggtitle("Product 3: endec_mu = 0.08, eff_len = 26")

mod3_plot_EIR <- df_mod3 %>%
  filter(ref == index_mod3) %>%
  ggplot()+
  aes(x = t, y = EIR_tot)+
  geom_line()+
  ylim(0, 200)+
  ggtitle("Product 3: endec_mu = 0.08, eff_len = 26")

mod3_plot_natural_dead_endec <- df_mod3 %>%
  filter(ref == index_mod3) %>%
  ggplot()+
  aes(x = t, y = mvx_dead_natural)+
  geom_line()+
  ylim(0, 200)+
  ggtitle("Product 3: endec_mu = 0.08, eff_len = 26")


comparison_products <- cowplot::plot_grid(mod1_mosq, mod2_mosq,mod3_plot_mosq,
                   mod1_killed_mosq, mod2_killed_mosq,mod3_plot_mosq_killed,
                   mod1_EIR, mod2_EIR, mod3_plot_EIR,
                   ncol =3,
                   nrow =3)




ggsave(comparison_products, file = "analysis/target-profiles-distrib-strat/comparison_products.svg")



cowplot::plot_grid(mod1_killed_mosq, mod2_killed_mosq,
                   mod3_plot_mosq_killed, nrow = 3)

#fitting
eff_len_endec4 <- 40
df_mod4 <-  readRDS("analysis/target-profiles-distrib-strat/killing_mose_endec_mod4.rds")
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
  filter(ref == index_mod4)
max(best_fit_mod4$mvx_dead)#54
max(best_fit_mod3$mvx_dead)#60
max(best_fit_mod2$mvx_dead) #58
max(mod1_df$mvx_dead) #58
#rbind them together and plot on top of each other.

mod_1_compar <- mod1_df %>%
  mutate(product = "product A") %>%
  select(t, endec_mu, eff_len, product, mv, mvx_dead, EIR_tot, clin_inc0to5, slide_prev0to5, betaa, beta_larval, KL,mvx_dead_natural)

ggplot(mod_1_compar, aes(x = t, y = mvx_dead))+
  geom_line()

best_fit_mod2 <- best_fit_mod2 %>%
  mutate(product = "product B") %>%
  select(t, endec_mu, eff_len, product, mv, mvx_dead, EIR_tot, clin_inc0to5, slide_prev0to5, betaa, beta_larval, KL,mvx_dead_natural)


best_fit_mod3 <- best_fit_mod3 %>%
  mutate(product = "product C") %>%
  select(t, endec_mu, eff_len, product, mv, mvx_dead, EIR_tot, clin_inc0to5, slide_prev0to5, betaa, beta_larval, KL ,mvx_dead_natural)

best_fit_mod4 <- best_fit_mod4 %>%
  mutate(product = "product D") %>%
  select(t, endec_mu, eff_len, product, mv, mvx_dead, EIR_tot, clin_inc0to5, slide_prev0to5, betaa, beta_larval, KL, mvx_dead_natural)

all_products <- do.call("rbind", list(mod_1_compar, best_fit_mod2, best_fit_mod3, best_fit_mod4))

plot_endec_dead_mosq <- ggplot(all_products, aes(x = t, y = mvx_dead, col = as.factor(product)))+
  geom_line()+
  geom_vline(aes(xintercept = IVM_start[1]), col = "black")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec1), col = "red", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec2), col = "green", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec3), col = "blue", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec4), col = "purple", linetype = "dashed")+
  xlim(350, 600)

plot_natural_killed_mosq <- ggplot(all_products, aes(x = t, y = mvx_dead_natural, col = as.factor(product)))+
  geom_line()+
  geom_vline(aes(xintercept = IVM_start[1]), col = "black")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec1), col = "red", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec2), col = "green", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec3), col = "blue", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec4), col = "purple", linetype = "dashed")+
  xlim(350, 600)

plot_eir_tot <- ggplot(all_products, aes(x = t, y = EIR_tot, col = as.factor(product)))+
  geom_line()+
  geom_vline(aes(xintercept = IVM_start[1]), col = "black")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec1), col = "red", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec2), col = "green", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec3), col = "blue", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec4), col = "purple", linetype = "dashed")+
  xlim(350, 600)

plot_mv <- ggplot(all_products, aes(x = t, y = mv, col = as.factor(product)))+
  geom_line()+
  geom_vline(aes(xintercept = IVM_start[1]), col = "black")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec1), col = "red", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec2), col = "green", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec3), col = "blue", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec4), col = "purple", linetype = "dashed")+
  xlim(350, 600)

plot_betaa <- ggplot(all_products, aes(x = t, y = betaa, col = as.factor(product)))+
  geom_line()+
  geom_vline(aes(xintercept = IVM_start[1]), col = "black")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec1), col = "red", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec2), col = "green", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec3), col = "blue", linetype = "dashed")+
  geom_vline(aes(xintercept = IVM_start[1] + eff_len_endec4), col = "purple", linetype = "dashed")+
  xlim(350, 600)


comparisons_TC <- cowplot::plot_grid(plot_endec_dead_mosq, plot_natural_killed_mosq,plot_eir_tot,
                   plot_mv, plot_betaa, ncol = 2, nrow = 3 )

#efficacy
mod0_times_prod1 <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec1))

mod0_times_prod2 <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec2))

mod0_times_prod3 <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec3))

mod0_times_prod4 <- mod0_df %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4))

mod1_times <- mod_1_compar %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec1))

mod2_times <- best_fit_mod2 %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec2))

mod3_times <- best_fit_mod3 %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec3))

mod4_times <- best_fit_mod4 %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4))


##efficacy calc##
#take mean in given period and compare

eff_period <- 30*4

mod0_eff <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_period)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))


mod1_eff <- mod_1_compar%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_period)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod2_eff <- best_fit_mod2%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_period)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod3_eff <- best_fit_mod3%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_period)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod4_eff <- best_fit_mod4%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_period)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

#eff EIR###
(mod0_eff$mean_EIR-mod1_eff$mean_EIR)/mod0_eff$mean_EIR #0.228..
(mod0_eff$mean_EIR-mod2_eff$mean_EIR)/mod0_eff$mean_EIR #0.227...
(mod0_eff$mean_EIR-mod3_eff$mean_EIR)/mod0_eff$mean_EIR #0.235...
(mod0_eff$mean_EIR-mod4_eff$mean_EIR)/mod0_eff$mean_EIR #0.221...

#eff prev###
(mod0_eff$mean_prev-mod1_eff$mean_prev)/mod0_eff$mean_prev #0.1712..
(mod0_eff$mean_prev-mod2_eff$mean_prev)/mod0_eff$mean_prev #0.1716..
(mod0_eff$mean_prev-mod3_eff$mean_prev)/mod0_eff$mean_prev #0.1722..
(mod0_eff$mean_prev-mod4_eff$mean_prev)/mod0_eff$mean_prev #0.1675..

#eff mv###
(mod0_eff$mean_mv-mod1_eff$mean_mv)/mod0_eff$mean_mv #0.0619..
(mod0_eff$mean_mv-mod2_eff$mean_mv)/mod0_eff$mean_mv #0.0622..
(mod0_eff$mean_mv-mod3_eff$mean_mv)/mod0_eff$mean_mv #0.0647..
(mod0_eff$mean_mv-mod4_eff$mean_mv)/mod0_eff$mean_mv #0.0562..


#efficacy for the duration of the longest product####
mod0_eff_long <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))


mod1_eff_long <- mod_1_compar%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod2_eff_long <- best_fit_mod2%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod3_eff_long <- best_fit_mod3%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod4_eff_long <- best_fit_mod4%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

#efficacy in duration of the longest product#
(mod0_eff_long$mean_EIR-mod1_eff_long$mean_EIR)/mod0_eff_long$mean_EIR #0.5595179
(mod0_eff_long$mean_EIR-mod2_eff_long$mean_EIR)/mod0_eff_long$mean_EIR #0.5716691
(mod0_eff_long$mean_EIR-mod3_eff_long$mean_EIR)/mod0_eff_long$mean_EIR #0.5602652
(mod0_eff_long$mean_EIR-mod4_eff_long$mean_EIR)/mod0_eff_long$mean_EIR #0.4235812

#eff prev###
(mod0_eff_long$mean_prev-mod1_eff_long$mean_prev)/mod0_eff_long$mean_prev #0.1590387
(mod0_eff_long$mean_prev-mod2_eff_long$mean_prev)/mod0_eff_long$mean_prev #0.1618007
(mod0_eff_long$mean_prev-mod3_eff_long$mean_prev)/mod0_eff_long$mean_prev #0.1722..
(mod0_eff_long$mean_prev-mod4_eff_long$mean_prev)/mod0_eff_long$mean_prev #0.1573832

#eff mv###
(mod0_eff_long$mean_mv-mod1_eff_long$mean_mv)/mod0_eff_long$mean_mv #0.2229541
(mod0_eff_long$mean_mv-mod2_eff_long$mean_mv)/mod0_eff_long$mean_mv #0.2253552
(mod0_eff_long$mean_mv-mod3_eff_long$mean_mv)/mod0_eff_long$mean_mv #0.2288135
(mod0_eff_long$mean_mv-mod4_eff_long$mean_mv)/mod0_eff_long$mean_mv #0.1766791

#

mod0_times_prod1  <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec1)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod0_times_prod2 <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec2)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod0_times_prod3 <- mod0_df%>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec3)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod0_times_prod4 <- mod0_df %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod1_times <- mod_1_compar %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec1)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod2_times <- best_fit_mod2 %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec2)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod3_times <- best_fit_mod3 %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec3)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

mod4_times <- best_fit_mod4 %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec4)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            mean_inc = mean(clin_inc0to5),
            mean_mv = mean(mv))

#efficacy for the duration of the product
(mod0_times_prod1$mean_EIR - mod1_times$mean_EIR)/mod0_times_prod1$mean_EIR #0.5531054
(mod0_times_prod2$mean_EIR - mod2_times$mean_EIR)/mod0_times_prod2$mean_EIR #0.5766688
(mod0_times_prod3$mean_EIR - mod3_times$mean_EIR)/mod0_times_prod3$mean_EIR #0.5470573
(mod0_times_prod4$mean_EIR - mod4_times$mean_EIR)/mod0_times_prod4$mean_EIR #0.4235812

(mod0_times_prod1$mean_prev - mod1_times$mean_prev)/mod0_times_prod1$mean_prev #0.1400387
(mod0_times_prod2$mean_prev - mod2_times$mean_prev)/mod0_times_prod2$mean_prev #0.1383426
(mod0_times_prod3$mean_prev - mod3_times$mean_prev)/mod0_times_prod3$mean_prev #0.142075
(mod0_times_prod4$mean_prev - mod4_times$mean_prev)/mod0_times_prod4$mean_prev #0.1482846

(mod0_times_prod1$mean_inc - mod1_times$mean_inc)/mod0_times_prod1$mean_inc #0.2149413
(mod0_times_prod2$mean_inc - mod2_times$mean_inc)/mod0_times_prod2$mean_inc #0.1789803
(mod0_times_prod3$mean_inc - mod3_times$mean_inc)/mod0_times_prod3$mean_inc #0.2497882
(mod0_times_prod4$mean_inc - mod4_times$mean_inc)/mod0_times_prod4$mean_inc #0.2840877

(mod0_times_prod1$mean_mv - mod1_times$mean_mv)/mod0_times_prod1$mean_mv #0.3025108
(mod0_times_prod2$mean_mv - mod2_times$mean_mv)/mod0_times_prod2$mean_mv #0.3362034
(mod0_times_prod3$mean_mv - mod3_times$mean_mv)/mod0_times_prod3$mean_mv #0.2858121
(mod0_times_prod4$mean_mv - mod4_times$mean_mv)/mod0_times_prod4$mean_mv #0.1766791

#not really a detectable difference in reduction in EIR between the scenarios.
#few more mosquitoes are able to transmit but not many.


#rbinding
df_mod3_product <- df_mod3 %>%
  #filter(between(t, IVM_start[1]-10, IVM_start[1]+eff_len_endec3+90)) %>%
  mutate(product = "product C") %>%
  select(product, eff_len, t,mv, mvx_dead, EIR_tot, clin_inc0to5, slide_prev0to5, betaa, beta_larval, KL,Ivtot)


df_mod2_product <- df_mod2 %>%
  #filter(between(t, IVM_start[1]-10, IVM_start[1]+eff_len_endec3+90)) %>%
  mutate(product = "product B") %>%
  select(product, eff_len, t,mv, mvx_dead, EIR_tot, clin_inc0to5, slide_prev0to5,betaa, beta_larval, KL,Ivtot)


df_mod1_product <- mod1_df %>%
  #filter(between(t, IVM_start[1]-10, IVM_start[1]+eff_len_endec3+90)) %>%
  mutate(product = "product A", eff_len = 23) %>%
  select(product, eff_len, t,mv, mvx_dead, EIR_tot, clin_inc0to5, slide_prev0to5,betaa, beta_larval, KL,Ivtot)


products <- do.call("rbind", list(df_mod1_product,  df_mod2_product,df_mod3_product))
class(products)

mvx_dead_plot <- ggplot(products, aes(x = t, y = mvx_dead, col = as.factor(product)))+
  geom_line()+
  xlim(400, 700)+
  theme(legend.position = c(0.8, 0.3))
Ivtot_plot <-  ggplot(products, aes(x = t, y = Ivtot, col = as.factor(product)))+
  geom_line()+
  xlim(400, 700)+
  theme(legend.position = "none")
mv_plot <- ggplot(products, aes(x = t, y = mv, col = as.factor(product)))+
  geom_line()+
  xlim(400, 700)+
  theme(legend.position ="none")
eir_tot_plot <- ggplot(products, aes(x = t, y = EIR_tot, col = as.factor(product)))+
  geom_line()+
  theme(legend.position ="none")
inc_plot <- ggplot(products, aes(x = t, y = clin_inc0to5, col = as.factor(product)))+
  geom_line()+
  ylim(0,0.01)+
  xlim(400, 700)+
  theme(legend.position ="none")
prev_plot <- ggplot(products, aes(x = t, y = slide_prev0to5, col = as.factor(product)))+
  geom_line()+
  ylim(0.25,0.75)+
  xlim(400, 700)+
  theme(legend.position ="none")
betaa_plot <- ggplot(products, aes(x = t, y = betaa, col = as.factor(product)))+
  geom_line()+
  xlim(400, 700)+
  theme(legend.position ="none")
beta_larval_plot <- ggplot(products, aes(x = t, y = beta_larval, col = as.factor(product)))+
  geom_line()+
  xlim(400, 700)+
  theme(legend.position ="none")
KL_plot <- ggplot(products, aes(x = t, y = KL, col = as.factor(product)))+
  geom_line()+
  xlim(400, 700)+
  theme(legend.position ="none")


cowplot::plot_grid(mvx_dead_plot, mv_plot, Ivtot_plot, eir_tot_plot, inc_plot, prev_plot,
                   betaa_plot, beta_larval_plot, KL_plot,
                   nrow = 5, ncol = 2)

ggplot(mod0_df, aes(x = t, y = clin_inc0to5))+
  geom_line()+
  ylim(0,0.05)


#then do the averted calcs.

mod0_summary <- mod0_df %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec1)) %>%
  summarise(mean_EIR_tot = mean(EIR_tot),
            mean_inc05 = mean(clin_inc0to5),
            mean_prev = mean(slide_prev0to5),
            mean_mv = mean(mv))

mod1_summary <- df_mod1_product %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec1)) %>%
  summarise(mean_EIR_tot = mean(EIR_tot),
            mean_inc05 = mean(clin_inc0to5),
            mean_prev = mean(slide_prev0to5),
            mean_mv = mean(mv))
mod2_summary <- df_mod2_product %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec2)) %>%
  summarise(mean_EIR_tot = mean(EIR_tot),
            mean_inc05 = mean(clin_inc0to5),
            mean_prev = mean(slide_prev0to5),
            mean_mv = mean(mv))
mod3_summary <- df_mod3_product %>%
  filter(between(t, IVM_start[1], IVM_start[1]+eff_len_endec2)) %>%
  summarise(mean_EIR_tot = mean(EIR_tot),
            mean_inc05 = mean(clin_inc0to5),
            mean_prev = mean(slide_prev0to5),
            mean_mv = mean(mv))


(mod0_summary$mean_EIR_tot-mod1_summary$mean_EIR_tot)/mod0_summary$mean_EIR_tot #55% decrease
(mod0_summary$mean_EIR_tot-mod2_summary$mean_EIR_tot)/mod0_summary$mean_EIR_tot #88%
(mod0_summary$mean_EIR_tot-mod3_summary$mean_EIR_tot)/mod0_summary$mean_EIR_tot #9.0%

(mod0_summary$mean_inc05-mod1_summary$mean_inc05)/mod0_summary$mean_inc05 #21% decrease
(mod0_summary$mean_inc05-mod2_summary$mean_inc05)/mod0_summary$mean_inc05 #34%
(mod0_summary$mean_inc05-mod3_summary$mean_inc05)/mod0_summary$mean_inc05 #6.1%

(mod0_summary$mean_prev-mod1_summary$mean_prev)/mod0_summary$mean_prev #14% decrease
(mod0_summary$mean_prev-mod2_summary$mean_prev)/mod0_summary$mean_prev #14%
(mod0_summary$mean_prev-mod3_summary$mean_prev)/mod0_summary$mean_prev #13%

(mod0_summary$mean_mv-mod1_summary$mean_mv)/mod0_summary$mean_mv #30% decrease
(mod0_summary$mean_mv-mod2_summary$mean_mv)/mod0_summary$mean_mv #71%
(mod0_summary$mean_mv-mod3_summary$mean_mv)/mod0_summary$mean_mv #2.3%


df_mod4 <- readRDS(file = "analysis/target-profiles-distrib-strat/killing_mose_endec_mod4.rds")
ggplot(df_mod4, aes(x = t, y = mvx_dead))+
  geom_line()
#fitting
df_mod4_dead <- df_mod4 %>%
  filter(t == IVM_start[1]+eff_len_endec4) %>%
  select(t, mv, mvx_dead, ref) %>%
  group_by(ref) %>%
  summarise(tot_mvx_dead = mvx_dead, mv = mv)
range(df_mod4_dead$tot_mvx_dead)
mod4_list <- split(df_mod4_dead, f = df_mod4_dead$ref)
error <- numeric()
for (i in 1:length(endec_mu_in)){
  error <- c(error, sum(mod1_endec_killed$mvx_dead - mod4_list[[i]]$tot_mvx_dead)^2)
}
index <- which.min(error) #1
endec_mu_in[index] #gives 0
df_mod4_dead[1,]



#if an endectocide is longer lasting than the 23-day killing product with endec_mu = 0.01:
#it can afford to have a much lower endec_mu, to kill approximately the same number of mosquitoes over the new product's killing period
#shorter lasting endectocides would need to have a much higher killing

df_mod1_dat <- mod1_df %>%
  select(t, mv, mvx_dead, EIR_tot, slide_prev0to5, endec_mu) %>%
  mutate(eff_len =23)

df_mod2_dat <- df_mod2 %>%
  filter(ref == 14) %>%
  select(t, mv, mvx_dead, EIR_tot, slide_prev0to5, endec_mu)%>%
  mutate(eff_len = 7)

df_mod3_dat <- df_mod3 %>%
  filter(ref == 11) %>%
  select(t, mv, mvx_dead, EIR_tot, slide_prev0to5, endec_mu)%>%
  mutate(eff_len = 14)

df_mod4_dat <- df_mod4 %>%
  filter(ref == 10) %>%
  select(t, mv, mvx_dead, EIR_tot, slide_prev0to5, endec_mu) %>%
  mutate(eff_len = 180)

df_list <- list(df_mod1_dat, df_mod2_dat, df_mod3_dat, df_mod4_dat)
df_bind <- do.call(rbind, df_list)


mv_plot <- ggplot(df_bind, aes(x = t, y = mv))+
  geom_line()+
  xlim(2000, 3650)+
  facet_grid(~eff_len, labeller = label_both)

mvx_dead_plot <- ggplot(df_bind, aes(x = t, y = mvx_dead))+
  geom_line()+
  xlim(2000, 3650)+
  facet_grid(~eff_len, labeller = label_both)

eir_plot <- ggplot(df_bind, aes(x = t, y = EIR_tot))+
  geom_line()+
  xlim(2000, 3650)+
  facet_grid(~eff_len, labeller = label_both)

prev_plot <- ggplot(df_bind, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  xlim(2000, 3650)+
  facet_grid(~eff_len, labeller = label_both)

cowplot::plot_grid(mv_plot, mvx_dead_plot, eir_plot, prev_plot)

#some diagnostic plots
mod1_plot <- ggplot(mod1_df, aes(x = t, y = mvx_dead))+
  geom_line()+
  ggtitle("endec_mu = 0.09, eff_len = 23")

mod2_plot <- df_mod2 %>%
  #filter(ref == 14) %>%
  filter(ref %in% c(1, 60, 80,100)) %>%
  ggplot()+
  aes(x = t, y = mvx_dead, col = as.factor(endec_mu))+
  geom_line()+
  ggtitle("eff_len = 7")+
  guides(color = "none")

mod3_plot <- df_mod3 %>%
  #filter(ref == 14) %>%
  filter(ref %in% c(1, 60, 80,100)) %>%
  ggplot()+
  aes(x = t, y = mvx_dead, col = as.factor(endec_mu))+
  geom_line()+
  ggtitle("eff_len = 14")+
  guides(color = "none")

mod4_plot <- df_mod4 %>%
  #filter(ref == 14) %>%
  filter(ref %in% c(1, 60, 80,100)) %>%
  ggplot()+
  aes(x = t, y = mvx_dead, col = as.factor(endec_mu))+
  geom_line()+
  ggtitle("eff_len = 30")+
  guides(color = "none")

mod5_plot <- df_mod5 %>%
  #filter(ref == 14) %>%
  filter(ref %in% c(1, 60, 80,100)) %>%
  ggplot()+
  aes(x = t, y = mvx_dead, col = as.factor(endec_mu))+
  geom_line()+
  ggtitle("eff_len = 180")+
  theme(legend.position = "top")

cowplot::plot_grid(mod1_plot, mod2_plot, mod3_plot, mod4_plot, mod5_plot, ncol = 2, nrow = 3)

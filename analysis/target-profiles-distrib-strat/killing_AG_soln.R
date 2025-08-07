#killing mosquitoes with constant mosquito emergence rate

#use odin_model_malaria_exp_decay.R and set the waning parameter to 0 - no waning
#script to compare different model types
devtools::load_all()
require(tidyverse)
#how much does mu_h change across different values of Q0 and ivm_cov?

# Provide a value of the annual EIR for this model run
#init_EIR_vec <- c(2, 25, 100) #low - 2, moderate - 15, high - 120 --> Ellie: low = 2, med = 25, high = 100

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years, turn ivermectin on when nets are 6m, 1y, 2.5yo
time_period <- 365
mda_int <- 30

#ivm on when nets are 6 months old
#net_seq <- seq(100, 3650, by = 3*365) #if other interventions e.g. nets are on, will need to model IVM on time in relation to this
IVM_begin1 <- 30
IVM_start <- c(IVM_begin1) #just modelling one pulse

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

#itn_on <- 100 #introduce nets 100 days into simulation

#net_seq <- seq(365, 3650, by = 3*365)
runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}


eff_lenA <- 23
#here the hazard profile is 1 everyday
ivm_parmsA <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_lenA),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)

#fits to get endec_mu and wane
endec_mu_vec <-0
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_varA <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_A <- list()
for (i in seq_len(nrow(df_varA))){
  my_list_A[[i]] <- as.numeric(df_varA[i,])
}

#endectocide model
modA_run <-  function(data_in){
  endec_mu_in <- data_in[1]
  wane_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay3.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 365,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 800,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms1$ttt,
    eff_len = ivm_parms1$eff_len,
    haz = ivm_parms1$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parmsA$ivm_min_age,
    ivm_max_age = ivm_parmsA$ivm_max_age,
    IVRM_start = ivm_parmsA$IVRM_start,
    Q0 = 1,
    endec_mu = endec_mu_in,
    wane = wane_in,
    #endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_modA <- function(){
  modA_out_list <- lapply(my_list_A, modA_run)
  res_modA_out <- lapply(modA_out_list, runfun)
  modA_df <- do.call(rbind, sapply(1:(nrow(df_varA)), function(x){
    df <- as.data.frame(res_modA_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(modA_df)
}

df_modA <- my_sim_modA()

#define constant inputs
mu <- df_modA$mu[1] #this ia baseline mortality rate of mosquitos
M0 <- df_modA$mv[1] # initial mosquito density (at eqm)

ggplot(df_modA, aes(x = t, y= mv))+
  geom_line()+
  ylim(0,200)

delta_D <- 100 #define the target number of mosquitoes to kill by intervention
time_range <- c(0, 365)
delta_t_vec <- c(23, 20, 26, 40) #vector of time durations for killing

# --- BASELINE PSI CALCULATION ---
# This calculates the value of psi required for equal birth and death rates.
# Note, to account for limits in machine precision, a lower bound is placed on
# 1 - delta_D_target/M0 = epsilon, where epsilon is the smallest non-zero
# normalised floating-point number. Also applied when delta_D_target >= M0.
baseline_psi <- function(Delta_t, delta_D_target, M0) {
  log_value <- max((M0 - delta_D_target) / M0, .Machine$double.xmin)
  -log(log_value) / Delta_t # Return value of psi
}

# --- FITTING FUNCTION FOR PSI ---
# This section calculates psi for a given a constant emergence rate.

# This is the objective function to minimise
psi_objective_function <- function(psi, Delta_t, delta_D_target, mu_v, M0) {
  term1 <- mu_v * Delta_t
  term2 <- (psi / (psi + mu_v)) * (1 - exp(-(psi + mu_v) * Delta_t))
  D_pred <- (psi * M0 / (psi + mu_v)) * (term1 + term2)
  (delta_D_target - D_pred)^2
}

# This fits the value of psi that minimises the objective function when constant emergence is true
fit_psi <- function(Delta_t, delta_D_target, mu_v, M0) {

  # Use the value of gamma under equal birth and death rates as initial guess
  guess <- baseline_psi(Delta_t, delta_D_target, M0)

  # Minimise the objective function
  fit <- optim(
    par = guess,
    fn = psi_objective_function,
    method = "L-BFGS-B",
    lower = 1e-6,
    upper = 100,
    control = list(factr = 1e2, pgtol = 1e-8, maxit = 1000),
    Delta_t = Delta_t,
    delta_D_target = delta_D_target,
    mu_v = mu_v,
    M0 = M0
  )

  # Return the estimate of psi
  fit$par

}

#endec1_int_rate = 0.09

gamma_vec <- numeric(length(delta_t_vec))
names(gamma_vec) <- paste0("Delta_t_", delta_t_vec)

for (i in seq_along(delta_t_vec)) {
  Delta_t <- delta_t_vec[i]
  gamma <- fit_gamma(Delta_t, delta_D, mu, M0)
  gamma_vec[i] <- gamma
}
gamma_vec
#this gives out gamma vecs which you can put into the subsequent models

#so then for the different endectocides...
#fits to get endec_mu and wane
endec_mu_vec <-gamma_vec[1]
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var1 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_1 <- list()
for (i in seq_len(nrow(df_var1))){
  my_list_1[[i]] <- as.numeric(df_var1[i,])
}

eff_len1 <- delta_t_vec[1]
#here the hazard profile is 1 everyday
ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len1),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)


#endectocide model
mod1_run <-  function(data_in){
  endec_mu_in <- data_in[1]
  wane_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay3.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 365,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 800,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms1$ttt,
    eff_len = ivm_parms1$eff_len,
    haz = ivm_parms1$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms1$ivm_min_age,
    ivm_max_age = ivm_parms1$ivm_max_age,
    IVRM_start = ivm_parms1$IVRM_start,
    Q0 = 1,
    endec_mu = endec_mu_in,
    wane = wane_in,
    #endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_mod1 <- function(){
  mod1_out_list <- lapply(my_list_1, mod1_run)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var1)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod1_df)
}

max(df_mod1$mv)

df_mod1 <- my_sim_mod1()
ggplot(df_mod1, aes(x = t, y  = mv))+
  geom_line()
ggplot(df_mod1, aes(x = t, y  = mvx_dead))+
  geom_line()


#next endec
#so then for the different endectocides...
#fits to get endec_mu and wane
endec_mu_vec <-gamma_vec[4]
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var2 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_2 <- list()
for (i in seq_len(nrow(df_var2))){
  my_list_2[[i]] <- as.numeric(df_var2[i,])
}



eff_len2 <- delta_t_vec[4]
#here the hazard profile is 1 everyday
ivm_parms2 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len2),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)


#endectocide model
mod2_run <-  function(data_in){
  endec_mu_in <- data_in[1]
  wane_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay3.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 365,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 300,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms1$ttt,
    eff_len = ivm_parms1$eff_len,
    haz = ivm_parms1$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms2$ivm_min_age,
    ivm_max_age = ivm_parms2$ivm_max_age,
    IVRM_start = ivm_parms2$IVRM_start,
    Q0 = 1,
    endec_mu = endec_mu_in,
    wane = wane_in,
    #endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_mod2 <- function(){
  mod2_out_list <- lapply(my_list_2, mod2_run)
  res_mod2_out <- lapply(mod2_out_list, runfun)
  mod2_df <- do.call(rbind, sapply(1:(nrow(df_var2)), function(x){
    df <- as.data.frame(res_mod2_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod2_df)
}

df_mod2 <- my_sim_mod2()
ggplot(df_mod2, aes(x = t, y  = mv))+
  geom_line()
ggplot(df_mod2, aes(x = t, y  = mvx_dead))+
  geom_line()
max(df_mod2$mvx_dead)
max(df_mod1$mvx_dead)

ggplot() +
  geom_line(data = df_mod1, aes(x = t, y = mvx_dead), color = "black") +
  geom_line(data = df_mod2, aes(x = t, y = mvx_dead), color = "blue")+
  ylim(0,100)


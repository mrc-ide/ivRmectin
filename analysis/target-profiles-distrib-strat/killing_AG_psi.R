devtools::load_all()
library(ggplot2)
library(tidyr)
library(dplyr)

# --- OPTIONS ---
# Set TRUE/FALSE options here
constant_emergence <- TRUE # phi = mu*M0 when TRUE or mu*M when FALSE
plot_R0_equals_1 <- FALSE
plot_Re_equals_1 <- FALSE
show_plot <- TRUE
save_plot <- TRUE
append_timestamp <- TRUE



#calculcate the endemic equilibrium in the absence of interventions
time_period <- 365*3
init_EIR <- 300
eff_lenA <- 23
IVM_start <- 365
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
endec_mu_vec <-0 #no intervention
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_varA <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_A <- list()
for (i in seq_len(nrow(df_varA))){
  my_list_A[[i]] <- as.numeric(df_varA[i,])
}

int_seq <- seq(IVM_start, IVM_start+24, by = 1)

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE, rtol = 1e-10, atol = 1e-12)
  modx <- mod$run(t = 0:time_period, tcrit = int_seq)
  op<- mod$transform_variables(modx)
  return(op)
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
    init_EIR = init_EIR,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parmsA$ttt,
    eff_len = ivm_parmsA$eff_len,
    haz = ivm_parmsA$haz,
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
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, theta2, Svtot, FOIv, clin_inc0to5, D_kill, betaa))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(modA_df)
}

df_modA <- my_sim_modA()
df_modA$mv[1]
df_modA$D_kill


# --- MODEL PARAMETERS ---
# Define model parameters here
params_base <- list(
  #N0 = 1000,              # Initial human population (should remain constant)

  M0 = df_modA$mv[1],              # Initial mosquito population

  #mu_h = 1 / (70 * 365),  # Human natural death rate
  ## (WHO Life Tables; 70y life expectancy)
  #
  #sigma_h = 1 / 12,       # Human incubation rate
  ## (~12 days incubation; Smith et al. 2012)
  #
  #gamma_h = 1 / 60,       # Human recovery rate
  ## (~60 days infectious; Griffin et al. 2010)
  #
  #omega_h = 1 / 180,      # Human loss of immunity rate
  ## (~180 days immune; White et al. 2014)
  #
  #beta_hv = 0.3,          # Daily transmission prob: human to mosquito
  ## (Smith et al. 2012)
  #
  mu_v = df_modA$mu[1],        # Mosquito natural death rate
  ## (~10 days lifespan; Lines et al. 1987)
  #
  #sigma_v = 1 / 10,       # Mosquito incubation rate (EIP)
  ## (~10 days at 25°C; Guerra et al. 2010)
  #
  #beta_vh = 0.3,          # Daily transmission prob: mosquito to human
  ## (Smith et al. 2012)

  #tau = 100,              # Start time of intervention (days)
  delta_D = 50,          # Target mosquitoes killed
  #times = seq(0, 500, by = 0.1), # All time points inclusive of pre-intervention
  delta_t_vec = c(23, 26, 30, 40, 50) # Intervention durations (days)
)
#params_base$m0 <- with(params_base, {M0 / N0}) # Initial mosquito to human ratio

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




# This fits the value of psi that minimises the objective function
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



psi_vec <- numeric(length(params_base$delta_t_vec))  # store

for (i in seq_along(sort(params_base$delta_t_vec))) {
  Delta_t <- sort(params_base$delta_t_vec)[i]

  if (constant_emergence) {
    psi_fit <- fit_psi(Delta_t, params_base$delta_D, params_base$mu_v, params_base$M0)
  } else {
    psi_fit <- baseline_psi(Delta_t, params_base$delta_D, params_base$M0)
  }

  psi_vec[i] <- psi_fit  # store the result
}

round(psi_vec, 5)

fit_psi(10, 50, 0.132, 103.4085)
fit_psi(50, 50, 0.132, 103.4085)

#overwrite it
##input from elsewhere
#psi_vec <- c(0.058196096, 0.027494417, 0.010318825, 0.005005920, 0.002465388)
int_info <- data.frame(params_base$delta_t_vec, psi_vec)
names(int_info) <- c("dur", "psi")


#killing duration 1####


#then run for different eff_len and corresponding psi_vec
eff_len1 <-int_info$dur[1]
#here the hazard profile is 1 everyday
ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len1),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)

#fits to get endec_mu and wane
endec_mu_vec <-int_info$psi[1]
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var1 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_1 <- list()
for (i in seq_len(nrow(df_var1))){
  my_list_1[[i]] <- as.numeric(df_var1[i,])
}

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
    init_EIR = init_EIR,
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
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, Svtot, mvx_dead_natural, FOIv, Evtot, Ivtot, clin_inc0to5, D_kill, betaa))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod1_df)
}
####

#killing duration 2####


#then run for different eff_len and corresponding psi_vec
eff_len2 <-int_info$dur[2]
#here the hazard profile is 1 everyday
ivm_parms2 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len2),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)

#fits to get endec_mu and wane
endec_mu_vec <-int_info$psi[2]
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var2 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_2 <- list()
for (i in seq_len(nrow(df_var2))){
  my_list_2[[i]] <- as.numeric(df_var2[i,])
}

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
    init_EIR = init_EIR,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms2$ttt,
    eff_len = ivm_parms2$eff_len,
    haz = ivm_parms2$haz,
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
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, Svtot,  mvx_dead_natural, FOIv, Evtot, Ivtot, clin_inc0to5, D_kill, betaa))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod2_df)
}
####

#killing duration 3####


#then run for different eff_len and corresponding psi_vec
eff_len3 <-int_info$dur[3]
#here the hazard profile is 1 everyday
ivm_parms3 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len3),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)

#fits to get endec_mu and wane
endec_mu_vec <-int_info$psi[3]
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var3 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_3 <- list()
for (i in seq_len(nrow(df_var3))){
  my_list_3[[i]] <- as.numeric(df_var3[i,])
}

#endectocide model
mod3_run <-  function(data_in){
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
    init_EIR = init_EIR,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms1$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start,
    Q0 = 1,
    endec_mu = endec_mu_in,
    wane = wane_in,
    #endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_mod3 <- function(){
  mod3_out_list <- lapply(my_list_3, mod3_run)
  res_mod3_out <- lapply(mod3_out_list, runfun)
  mod3_df <- do.call(rbind, sapply(1:(nrow(df_var3)), function(x){
    df <- as.data.frame(res_mod3_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, Svtot, mvx_dead_natural, FOIv, Evtot, Ivtot, clin_inc0to5, D_kill, betaa))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod3_df)
}
####

#killing duration 4####


#then run for different eff_len and corresponding psi_vec
eff_len4 <-int_info$dur[4]
#here the hazard profile is 1 everyday
ivm_parms4 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len4),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)

#fits to get endec_mu and wane
endec_mu_vec <-int_info$psi[4]
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var4 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_4 <- list()
for (i in seq_len(nrow(df_var4))){
  my_list_4[[i]] <- as.numeric(df_var4[i,])
}

#endectocide model
mod4_run <-  function(data_in){
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
    init_EIR = init_EIR,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms4$ttt,
    eff_len = ivm_parms4$eff_len,
    haz = ivm_parms4$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms4$ivm_min_age,
    ivm_max_age = ivm_parms4$ivm_max_age,
    IVRM_start = ivm_parms4$IVRM_start,
    Q0 = 1,
    endec_mu = endec_mu_in,
    wane = wane_in,
    #endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_mod4 <- function(){
  mod4_out_list <- lapply(my_list_4, mod4_run)
  res_mod4_out <- lapply(mod4_out_list, runfun)
  mod4_df <- do.call(rbind, sapply(1:(nrow(df_var4)), function(x){
    df <- as.data.frame(res_mod4_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, Svtot, mvx_dead_natural, FOIv, Evtot, Ivtot, clin_inc0to5, D_kill, betaa))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod4_df)
}
####

#killing duration 5####


#then run for different eff_len and corresponding psi_vec
eff_len5 <-int_info$dur[5]
#here the hazard profile is 1 everyday
ivm_parms5 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len5),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)

#fits to get endec_mu and wane
endec_mu_vec <-int_info$psi[5]
wane <-  0
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var5 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_5 <- list()
for (i in seq_len(nrow(df_var5))){
  my_list_5[[i]] <- as.numeric(df_var5[i,])
}

#endectocide model
mod5_run <-  function(data_in){
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
    init_EIR = init_EIR,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms5$ttt,
    eff_len = ivm_parms5$eff_len,
    haz = ivm_parms5$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms5$ivm_min_age,
    ivm_max_age = ivm_parms5$ivm_max_age,
    IVRM_start = ivm_parms5$IVRM_start,
    Q0 = 1,
    endec_mu = endec_mu_in,
    wane = wane_in,
    #endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_mod5 <- function(){
  mod5_out_list <- lapply(my_list_5, mod5_run)
  res_mod5_out <- lapply(mod5_out_list, runfun)
  mod5_df <- do.call(rbind, sapply(1:(nrow(df_var5)), function(x){
    df <- as.data.frame(res_mod5_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, Svtot, mvx_dead_natural, FOIv, Evtot, Ivtot, clin_inc0to5, D_kill, betaa))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod5_df)
}
####


#run all the models
df_mod1 <- my_sim_mod1()
df_mod2 <- my_sim_mod2()
df_mod3 <- my_sim_mod3()
df_mod4 <- my_sim_mod4()
df_mod5 <- my_sim_mod5()

all_products <- do.call("rbind", list(df_mod1, df_mod2, df_mod3, df_mod4, df_mod5))

all_products %>%
  group_by(endec_mu) %>%
  summarise(max_dead = max(mvx_dead),
            prop_killed = max(mvx_dead)/mv[1]) #differing by max 9 mosquitoes. Check against Andrew's model



all_products_info <- left_join(all_products, int_info, by = c("endec_mu" = "psi"))

all_products_info <- all_products_info %>%
  mutate(killing_info = paste("duration =", dur, ", endec_mu =", round(endec_mu, 3)))

all_products_info %>%
  group_by(killing_info) %>%
  summarise(max_dead = max(mvx_dead),
            prop_killed = max(mvx_dead)/mv[1])

emergence_plot <- ggplot(all_products_info, aes(x = t, y = betaa, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  #geom_line( data = all_products_info,aes(x = t, y = mvx_dead, col = as.factor(killing_info)), linetype = "dashed")+
  labs(col = "Killing strategy")+
  theme(legend.position = c(0.7, 0.5))+
  #ylab("Number of intervention-killed \nmosquitoes")+
  guides(col = "none")+
  #xlim(300, 700)+
  geom_vline(xintercept = IVM_start, linetype = "dashed", col = "black")+
  geom_vline(xintercept = IVM_start+60, linetype = "dashed", col = "black")

ggplot(all_products_info, aes(x = t, y = D_kill, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  geom_line( data = all_products_info,aes(x = t, y = mvx_dead, col = as.factor(killing_info)), linetype = "dashed")+
  labs(col = "Killing strategy")+
  ylab("Number of intervention-killed \nmosquitoes")+
  #guides(col = "none")+
  #xlim(300, 700)+
  geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
  geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")

mosq_killed_plot <- ggplot(all_products_info, aes(x = t, y = mvx_dead, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  ylab("Number of intervention-killed \nmosquitoes")+
  guides(col = "none")+
  #xlim(300, 700)+
  geom_vline(xintercept = IVM_start, linetype = "dashed", col = "black")+
  geom_vline(xintercept = IVM_start+60, linetype = "dashed", col = "black")

prop_killed <- ggplot(all_products_info, aes(x = t, y = mvx_dead/mv[1], col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  ylab("Percentage ntervention-killed \nmosquitoes")+
  guides(col = "none")+
  #xlim(300, 700)+
  ylim(0,1)+
  geom_vline(xintercept = IVM_start, linetype = "dashed", col = "black")+
  geom_vline(xintercept = IVM_start+60, linetype = "dashed", col = "black")

#need to run longer to get back to equilibrium??
prev_plot <- ggplot(all_products_info, aes(x = t, y = slide_prev0to5, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  ylim(0.8,0.85)+
  guides(col = "none")+
  #xlim(300, 700)+
  geom_vline(xintercept = IVM_start, linetype = "dashed", col = "black")+
  geom_vline(xintercept = IVM_start+60, linetype = "dashed", col = "black")

mosq_density_plot <- ggplot(all_products_info, aes(x = t, y = mv, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  #xlim(300, 700)+
  guides(col = "none")+
  geom_vline(xintercept = IVM_start, linetype = "dashed", col = "black")+
  geom_vline(xintercept = IVM_start+60, linetype = "dashed", col = "black")

ggplot(all_products_info, aes(x = t, y = Svtot, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  xlim(300, 700)+
  theme(legend.position = c(0.7,0.45))+
  geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
  geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")



Svtot_plot <- ggplot(all_products_info, aes(x = t, y = Svtot, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  xlim(300, 700)+
  theme(legend.position = c(0.7,0.45))+
  guides(col = "none")+
  geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
  geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")

ggplot(all_products_info, aes(x = t, y = FOIv, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  xlim(300, 700)+
  theme(legend.position = c(0.7,0.45))+
  guides(col = "none")+
  geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
  geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")

ggplot(all_products_info, aes(x = t, y = Evtot, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  xlim(300, 700)+
  theme(legend.position = c(0.7,0.45))+
  guides(col = "none")+
  geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
  geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")

ggplot(all_products_info, aes(x = t, y = Ivtot, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  xlim(300, 700)+
  theme(legend.position = c(0.7,0.45))+
  guides(col = "none")+
  geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
  geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")

ggplot(all_products_info, aes(x = t, y = mvx_dead_natural, col = as.factor(killing_info)))+
    geom_line(size = 1)+
    labs(col = "Killing strategy")+
    xlim(300, 700)+
    theme(legend.position = c(0.7,0.45))+
    #guides(col = "none")+
    geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
    geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")

FOIv_plot <- ggplot(all_products_info, aes(x = t, y = FOIv, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  #xlim(300, 700)+
  theme(legend.position = c(0.7,0.45))+
  ylim(0.00928, 0.0094)+
  guides(col = "none")+
  geom_vline(xintercept = IVM_start, linetype = "dashed", col = "black")+
  geom_vline(xintercept = IVM_start+60, linetype = "dashed", col = "black")

inc_plot <- ggplot(all_products_info, aes(x = t, y = clin_inc0to5, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  #xlim(300, 700)+
  theme(legend.position = c(0.7,0.7))+
  ylim(0, 0.01)+
  guides(col = "none")+
  geom_vline(xintercept = IVM_start, linetype = "dashed", col = "black")+
  geom_vline(xintercept = IVM_start+60, linetype = "dashed", col = "black")

cowplot::plot_grid(mosq_killed_plot, prev_plot, mosq_density_plot,inc_plot,
                   nrow = 4)

cowplot::plot_grid(emergence_plot, mosq_density_plot,
                   mosq_killed_plot, prop_killed,
                   FOIv_plot, prev_plot, inc_plot)


#all killing similar proportion of mosquitoes
all_products %>%
  group_by(endec_mu) %>%
  summarise(max_dead = max(mvx_dead),
            prop_killed = max(mvx_dead)/mv[1])


#tot cases in baseline
cases_baseline <- df_modA %>%
  filter(between(t, IVM_start, IVM_start+60)) %>%
  summarise(cases_per_1000 = sum(clin_inc0to5)*1000) #this is the total number in that 60 day period

#cases averted
impact_measure <- all_products_info %>%
  filter(between(t, IVM_start, IVM_start+60)) %>%
  group_by(killing_info)%>%
  summarise(max_dead = max(mvx_dead),
            prop_killed = max(mvx_dead)/mv[1],
            cases_averted = sum(clin_inc0to5*1000)/cases_baseline$cases_per_1000)

#see a non-linear relationship between the total number of cases averted
ggplot(impact_measure, aes(x = prop_killed, y = cases_averted, col = as.factor(killing_info)))+
  geom_point()

#repeat this as a point estimate

cases_baseline_point <- df_modA %>%
  filter(t == IVM_start+60) %>%
  summarise(cases_per_1000 = clin_inc0to5*1000,
            prev = slide_prev0to5*100)

all_products_info %>%
  filter(t == IVM_start+60) %>%
  group_by(killing_info)%>%
  summarise(max_dead = mvx_dead,
            prop_killed = mvx_dead/mv[1],
            cases_averted_percent = ((cases_baseline_point$cases_per_1000 - (clin_inc0to5*1000))/cases_baseline_point$cases_per_1000)*100,
            cases_averted_abs = cases_baseline_point$cases_per_1000 - (clin_inc0to5*1000),
            inc = clin_inc0to5*1000,
            prev = slide_prev0to5*100)

#Carlos' output is incidence: see much bigger difference in incidence across the scenarios but very small difference in prevalence.
#I think the different interpretations come down to dynamics

results_long <- readRDS("C:/Users/nc1115/Documents/github/how_to_get_away_with_murder/results_long.rds")

D_AG <- results_long %>%
  filter(Compartment == "D") %>%
  rename(endec_mu = psi) %>%
  select(time, Value, endec_mu) %>%
  rename(D = Value) %>%
  mutate(model = "AG")

mv_AG <- results_long %>%
  filter(Compartment == "M") %>%
  rename(endec_mu = psi) %>%
  select(time, Value, endec_mu) %>%
  rename(M = Value) %>%
  mutate(model = "AG")



ggplot(all_products_info, aes(x = t, y = mv, col = as.factor(killing_info)))+
  geom_line(size = 1)+
  labs(col = "Killing strategy")+
  #xlim(300, 700)+
  theme(legend.position = c(0.7,0.45))+
  geom_vline(xintercept = 365, linetype = "dashed", col = "black")+
  geom_vline(xintercept = 365+60, linetype = "dashed", col = "black")

D_NC <-  all_products_info %>%
  select(t, mvx_dead, endec_mu) %>%
  rename(time = t,
         D = mvx_dead) %>%
  mutate(model = "NC")

mv_NC <- all_products_info %>%
  select(t, mv, endec_mu) %>%
  rename(time = t,
         M = mv) %>%
  mutate(model = "NC")

mv_models <- rbind(mv_AG, mv_NC) %>%
  mutate(endec_mu = round(endec_mu, 5))

ggplot(mv_models, aes(x = time, y = M, col = as.factor(endec_mu), linetype = as.factor(model)))+
  geom_line()

mv_models_wide <- mv_models %>%
  pivot_wider(names_from = model, values_from = M) %>%
  mutate(diff = AG-NC)

ggplot(mv_models_wide, aes(x = time, y = diff, col = as.factor(endec_mu)))+
  geom_line()

D_models <- rbind(D_AG, D_NC) %>%
  mutate(endec_mu = round(endec_mu, 5))

ggplot(D_models, aes(x = time, y = D, col = as.factor(endec_mu), linetype = as.factor(model)))+
  geom_line(size = 1.1)+
  theme_bw()

#number of killed mosquitoes seems to scale with difference in mosquito densities

D_models %>%
  group_by(model, endec_mu) %>%
  summarise(D_kill = max(D)) %>%
  pivot_wider(names_from = model, values_from = D_kill) %>%
  mutate(diff_int_killed = AG-NC)

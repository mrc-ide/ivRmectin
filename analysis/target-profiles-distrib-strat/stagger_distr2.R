#exploring staggered distribution

#model: odin_model_endectocide_staggered.R

devtools::load_all()
require(tidyverse)

time_period <- 365*2 #long run to get to eqm
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
start <- (365*2)+200 #time for starting distr, help get to reach eqm, and to distr at right time when seasonality introduced
num_distr <- 3 #either all in one, or over 3 distributions
target_cov <- 0.9
stag_cov <- target_cov/num_distr # this is what is feasible for field teams

#over 30 days
stag_30 <- 30
times_30_stag <- seq(start, start+ stag_30, length.out = num_distr)

#the start time for the different subpops
stag_30_IVM_start_1 <- times_30_stag[1]
stag_30_IVM_start_2 <- times_30_stag[2]
stag_30_IVM_start_3 <- times_30_stag[3]

#set IVM params for Hannah's mosq model with hazards#
ivm_parms_30_stag <- ivRmectin::ivm_fun_stag(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times_1 = stag_30_IVM_start_1,
  IVM_start_times_2 = stag_30_IVM_start_2,
  IVM_start_times_3 = stag_30_IVM_start_3,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

#over 10d
stag_10 <- 10
times_10_stag <- seq(start, start+ stag_10, length.out = num_distr)

#the start time for the different subpops
stag_10_IVM_start_1 <- times_10_stag[1]
stag_10_IVM_start_2 <- times_10_stag[2]
stag_10_IVM_start_3 <- times_10_stag[3]

#set IVM params for Hannah's mosq model with hazards#
ivm_parms_10_stag <- ivRmectin::ivm_fun_stag(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times_1 = stag_10_IVM_start_1,
  IVM_start_times_2 = stag_10_IVM_start_2,
  IVM_start_times_3 = stag_10_IVM_start_3,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

#QO for stephensi, arabiensis, gambiae and arabiensis-like vectors
Q0_in <- c(0.92)
init_EIR_in <- c(100)

df_var_stag <- expand.grid(Q0 = Q0_in, ivm_cov = stag_cov, init_EIR = init_EIR_in)

my_list_stag <- list()
for (i in seq_len(nrow(df_var_stag))){
  my_list_stag[[i]] <- as.numeric(df_var_stag[i,])
}

a <- Sys.time()

mod_30d <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_staggered2.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_30_stag$ttt,
    eff_len = ivm_parms_30_stag$eff_len,
    haz = ivm_parms_30_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_30_stag$ivm_min_age,
    ivm_max_age = ivm_parms_30_stag$ivm_max_age,
    IVRM_start_1 = ivm_parms_30_stag$IVRM_start_1,
    IVRM_start_2 = ivm_parms_30_stag$IVRM_start_2,
    IVRM_start_3 = ivm_parms_30_stag$IVRM_start_3,
    Q0 = Q0_in,
    num_distr = 3
  )
  return(output)
}

my_sim_mod_30d <- function(){
  mod_out_list <- lapply(my_list_stag, mod_30d)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_stag)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvtot_1, mvtot_2,mvtot_3, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "stag-30d"))}, simplify = F))
  return(mod_df)
}

df_mod30d <- my_sim_mod_30d()

ggplot(df_mod30d, aes(x = t, y = mv))+
  geom_line()

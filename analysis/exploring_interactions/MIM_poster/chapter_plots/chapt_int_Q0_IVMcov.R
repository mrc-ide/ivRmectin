devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
require(tidyverse)
bites_Bed_vector <- c(0.25, 0.5, 0.75, 0.95)
Q0_vector <-c(0.25, 0.5, 0.75, 0.95)
res_vector <- c(0, 0.1, 0.5, 0.7, 0.9)
itn_cov_vector <- c(0.2, 0.8, 0.2)
ivm_cov_vector <- c(0.3,0.5,0.7,0.9)

#time_period <- 3650 # run model for 10 years, turn ivermectin on when nets are 6m, 1y, 2.5yo
time_period <- 365*10
mda_int <- 30
ivm_cov = 0.7
#itn_cov_in = 0.8

itn_on <- 100 #introduce nets 100 days into simulation

net_seq <- seq(100, 3650, by = 3*365)

#ivm on when nets are 6 months old
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)


ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocide-specific parameters
source("R/mda_ivm_functions.R")

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#IVM_starters <- c(IVM_start1[1], IVM_start2[1], IVM_start3[1])

IVM_start <- numeric(length = 3)


#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivRmectin::ivm_fun(
  #IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start1,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

eff_len <- 23

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance %in% res_vector)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)


pyr_only_d_ITN0 <- df_pyr_only$dn0_med

pyr_param_list_og <- list()
pyr_param_list <- list()

Q0_ivm_grid <- expand.grid(bites_Bed = bites_Bed_vector[4], dn0_med = df_pyr_only$dn0_med[1], init_EIR = 100, itn_cov = 0.6, Q0 = Q0_vector,
                       ivm_cov = ivm_cov_vector)


Q0_ivm_df <- left_join(Q0_ivm_grid,
                       df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                       by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

Q0_ivm_list <- list()
for (i in seq_len(nrow(Q0_ivm_df))){
  Q0_ivm_list[[i]] <- as.numeric(Q0_ivm_df[i,])
}


antag_ITN_cov_loop <- function(itn_type_ivm_param){
  bites_Bed_in <- itn_type_ivm_param[1]
  d_ITN0_in <- itn_type_ivm_param[2]
  init_EIR_in <- itn_type_ivm_param[3]
  itn_cov_in <-itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[7]
  itn_half_life_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = itn_on,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}
my_sim_antag_ITN_bb_Q0 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN <- lapply(Q0_ivm_list, antag_ITN_cov_loop)
  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(Q0_ivm_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5, ivm_cov))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "bb-Q0", int = "ITN"))}, simplify = F))
  return(pyr_out_df_antag_ITN)

}

antag_ITN_bb_Q0 <- my_sim_antag_ITN_bb_Q0()
saveRDS(antag_ITN_bb_Q0, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/Q0_ITN_only.rds")


#antag with IVM
antag_ITN_IVM_cov_loop <- function(itn_type_ivm_param){
  bites_Bed_in <- itn_type_ivm_param[1]
  d_ITN0_in <- itn_type_ivm_param[2]
  init_EIR_in <- itn_type_ivm_param[3]
  itn_cov_in <-itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  ivm_cov_in <- itn_type_ivm_param[6]
  r_ITN0_in <- itn_type_ivm_param[7]
  itn_half_life_in <- itn_type_ivm_param[8]
  #IVRM_start_in <- ivm_nets_starting
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = itn_on,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_antag_ITN_IVM_bb_res <- function(){
  #pyr_out_list_antag_ITN_IVM <- purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN_IVM <- lapply(Q0_ivm_list, antag_ITN_IVM_cov_loop)
  res_pyr_out_antag_ITN_IVM <- lapply(pyr_out_list_antag_ITN_IVM, runfun) #put these values into the model
  pyr_out_df_antag_ITN_IVM <- do.call(rbind, sapply(1:(nrow(Q0_ivm_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN_IVM[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5, ivm_cov))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "Q0-ivm-cov", int = "ITN_IVM"))}, simplify = F))
  return(pyr_out_df_antag_ITN_IVM)

}

antag_ITN_IVM_bb_res <- my_sim_antag_ITN_IVM_bb_res()
saveRDS(antag_ITN_IVM_bb_res, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/Q0_ITN_IVM_cov.rds")

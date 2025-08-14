#bites Bed analysis

#model runs for the exploring interactions paper

#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
require(tidyverse)
bites_Bed_vector <- c(0.25, 0.5, 0.75, 0.95)
Q0_vector <-c(0.25, 0.5, 0.75, 0.95)
res_vector <- c(0, 0.1, 0.5, 0.7, 0.9)
itn_cov_vector <- seq(0.2, 0.8, 0.2)

#4 figures:
#1) Dynamics plot for A.gambiae-like vector at 10% resistance in different transmission settings
#2) Efficacy plot: predicted by anatagonisitic and additive model, shapes and colours for different phi-B, facet by Q0 and transmission setting
#3) Efficacy by species: creat a geom tile, for each transmission setting, and show the relative different in EIR and prevalence. Y is resistance, x is species

#relatives are LLIN & IVM compared to LLIN only
#ivRmectin model has 0.89 as bites_Bed gambiae, going to reset as 0.85 from PNAS paper

# Create a vector of age categories for the model
#init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR_vec <- c(2, 25, 100) #low - 2, moderate - 15, high - 120 --> Ellie: low = 2, med = 25, high = 100

# Provide the length of time (in days) that you want to run the model for
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

#when nets are 1yo

#when nets are 1yo

IVM_begin2 <- net_seq[3]+(1*365)
IVM_start2 <- c(IVM_begin2, IVM_begin2+mda_int, IVM_begin2 + mda_int + mda_int)

y2.5 <- (365*2.5)
#when nets are 2.5yo
IVM_begin3 <- net_seq[3]+y2.5
IVM_start3 <- c(IVM_begin3, IVM_begin3+mda_int, IVM_begin3+mda_int+mda_int)


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

IVM_starting <- list(IVM_start1, IVM_start2, IVM_start3)
#IVM_starting <- list(IVM_start1) #just running with 6m nets, where antagonism highest


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
#ivm_on <- IVM_start[1]
#ivm_off <- IVM_start[3]+eff_len



path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance == 0 | resistance  == 0.1 | resistance == 0.5 | resistance == 0.7 | resistance == 0.9)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

pyr_only_d_ITN0 <- df_pyr_only$dn0_med

pyr_param_list_og <- list()
pyr_param_list <- list()

#repeat for pyrethroid nets
#pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vector, bites_Bed = bites_Bed_vector,
 #                                init_EIR = init_EIR_vec, Q0 = Q0_vector)

bb_res_grid <- expand.grid(bites_Bed = bites_Bed_vector, dn0_med = pyr_only_d_ITN0, init_EIR = init_EIR_vec, itn_cov = 0.8, Q0 = Q0_vector[4])
bb_cov_grid <- expand.grid(bites_Bed =bites_Bed_vector,dn0_med = df_pyr_only$dn0_med[1],init_EIR = init_EIR_vec,itn_cov = itn_cov_vector,  Q0 = Q0_vector[4])
bb_Q0_grid <- expand.grid(bites_Bed = bites_Bed_vector, dn0_med = df_pyr_only$dn0_med[1], init_EIR = init_EIR_vec, itn_cov = 0.8, Q0 = Q0_vector)

bb_res_df <- left_join(bb_res_grid,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)


#filter to make it really simple
bb_res_df <- bb_res_df %>%
  filter(init_EIR == 100 & itn_cov == 0.8) %>% #high transmission setting
  select(-resistance)
bb_res_df
dim(bb_res_df) #20,7
bb_res_list <- list()

for (i in seq_len(nrow(bb_res_df))){
  bb_res_list[[i]] <- as.numeric(bb_res_df[i,])
}

##
bb_cov_df <-left_join(bb_res_grid,
                      df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                      by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365, Q0 = 0.95) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)


#filter to make it really simple
bb_cov_df <- bb_cov_df %>%
  filter(init_EIR == 100 & itn_cov == 0.8) %>% #high transmission setting
  select(-resistance)
dim(bb_cov_df) #20.7

bb_cov_list <- list()

for (i in seq_len(nrow(bb_res_df))){
  bb_cov_list[[i]] <- as.numeric(bb_cov_df[i,])
}
##
bb_Q0_df <- left_join(bb_Q0_grid,
                       df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                       by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365, Q0 = 0.95, itn_cov = 0.8) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)


#filter to make it really simple
bb_Q0_df <- bb_Q0_df %>%
  filter(init_EIR == 100 & itn_cov == 0.8) %>% #high transmission setting
  select(-resistance)
dim(bb_Q0_df)
bb_Q0_list <- list()

for (i in seq_len(nrow(bb_Q0_df))){
  bb_Q0_list[[i]] <- as.numeric(bb_Q0_df[i,])
}

antag_ITN_cov_loop <- function(itn_type_ivm_param){
  bites_Bed_in <- itn_type_ivm_param[1]
  d_ITN0_in <- itn_type_ivm_param[2]
  init_EIR_in <- itn_type_ivm_param[3]
  itn_cov_in <-itn_type_ivm_param[4]
  Q0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
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

my_sim_antag_ITN_bb_res <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN <- lapply(bb_res_list, antag_ITN_cov_loop)
  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(bb_res_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "bb-res"))}, simplify = F))
  return(pyr_out_df_antag_ITN)

}

antag_ITN_bb_res <- my_sim_antag_ITN_bb_res()

my_sim_antag_ITN_bb_cov <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN <- lapply(bb_cov_list, antag_ITN_cov_loop)
  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(bb_cov_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "bb-cov"))}, simplify = F))
  return(pyr_out_df_antag_ITN)

}

antag_ITN_bb_cov <- my_sim_antag_ITN_bb_cov()

my_sim_antag_ITN_bb_Q0 <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN <- lapply(bb_Q0_list, antag_ITN_cov_loop)
  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(bb_Q0_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "bb-Q0"))}, simplify = F))
  return(pyr_out_df_antag_ITN)

}

antag_ITN_bb_Q0 <- my_sim_antag_ITN_bb_Q0()


#antag with IVM
antag_ITN_IVM_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  Q0_in <- itn_type_ivm_param[7]
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
    ivm_cov_par = ivm_parms1$ivm_cov_par, # proportion of popuulation receiving the endectocide
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
  pyr_out_list_antag_ITN_IVM <- lapply(bb_res_list, antag_ITN_IVM_cov_loop)
  res_pyr_out_antag_ITN_IVM <- lapply(pyr_out_list_antag_ITN_IVM, runfun) #put these values into the model
  pyr_out_df_antag_ITN_IVM <- do.call(rbind, sapply(1:(nrow(bb_res_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN_IVM[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "bb_res"))}, simplify = F))
  return(pyr_out_df_antag_ITN_IVM)

}

antag_ITN_IVM_bb_res <- my_sim_antag_ITN_IVM_bb_res()

my_sim_antag_ITN_IVM_bb_cov <- function(){
  #pyr_out_list_antag_ITN_IVM <- purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN_IVM <- lapply(bb_cov_list, antag_ITN_IVM_cov_loop)
  res_pyr_out_antag_ITN_IVM <- lapply(pyr_out_list_antag_ITN_IVM, runfun) #put these values into the model
  pyr_out_df_antag_ITN_IVM <- do.call(rbind, sapply(1:(nrow(bb_cov_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN_IVM[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "bb_cov"))}, simplify = F))
  return(pyr_out_df_antag_ITN_IVM)

}

antag_ITN_IVM_bb_cov <- my_sim_antag_ITN_IVM_bb_cov()

my_sim_antag_ITN_IVM_bb_Q0 <- function(){
  #pyr_out_list_antag_ITN_IVM <- purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN_IVM <- lapply(bb_Q0_list, antag_ITN_IVM_cov_loop)
  res_pyr_out_antag_ITN_IVM <- lapply(pyr_out_list_antag_ITN_IVM, runfun) #put these values into the model
  pyr_out_df_antag_ITN_IVM <- do.call(rbind, sapply(1:(nrow(bb_Q0_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN_IVM[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "bb_Q0"))}, simplify = F))
  return(pyr_out_df_antag_ITN_IVM)

}

antag_ITN_IVM_bb_cov <- my_sim_antag_ITN_IVM_bb_Q0()


#bound antag models
antag_LLIN <- antag_ITN %>%
  mutate(model = "antag_LLIN")
antag_LLIN_IVM <- antag_ITN_IVM %>%
  mutate(model = "antag_LLIN_IVM")

#antag <- rbind(antag_LLIN, antag_LLIN_IVM)
#add <- rbind(add_LLIN, add_LLIN_IVM)
saveRDS(antag, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/antag_bb_combos.rds")
#saveRDS(add, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/add.rds")


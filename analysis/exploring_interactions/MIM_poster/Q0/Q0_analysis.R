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
Q0_vector <- c(0.25, 0.5, 0.75, 0.95)
bites_Bed_vector <-       rep(0.95, 4)


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
ivm_cov = 0.8
itn_cov_in = 0.8

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
# Sourcing the extra functions required to generate the endectocidepecific parameters
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

res_vector <- c(0, 0.1, 0.5, 0.7, 0.9)
itn_cov_vector <- seq(0.2, 0.8, 0.2)

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
pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vector, Q0 = Q0_vector,
                                 init_EIR = init_EIR_vec)
pyr_param_df <- left_join(pyr_param_df_crit,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365, bites_Bed = 0.95) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)
head(pyr_param_df)

#filter to make it really simple

pyr_param_df <- pyr_param_df %>%
  filter(resistance == 0.0 & init_EIR == 25 & itn_cov == 0.8) %>%
  select(-resistance)
head(pyr_param_df) #anyway filtering to no resistance and 80% coverage

for (i in seq_len(nrow(pyr_param_df))){
  pyr_param_list[[i]] <- as.numeric(pyr_param_df[i,])
}

antag_ITN_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  Q0_in <- itn_type_ivm_param[3]
  #IVRM_start_in <- itn_type_ivm_param[10] #failed because was getting made into numeric
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

my_sim_antag_ITN <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN <- lapply(pyr_param_list, antag_ITN_cov_loop)
  res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model
  pyr_out_df_antag_ITN <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
  return(pyr_out_df_antag_ITN)

}

antag_ITN <- my_sim_antag_ITN()


#antag with IVM
antag_ITN_IVM_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  Q0_in <- itn_type_ivm_param[3]
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
    ivm_cov_par = 0.8, # proportion of popuulation receiving the endectocide
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


#x <- rep(ivm_nets_starting, 2) #do reps to ensure have the same length
#y <- rep(pyr_param_list[1:2], 3)
##
#z <-purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
#z2 <- lapply(z, runfun) #put these values into the model
#z2df  <- do.call(rbind,
#                 sapply(1:(6), function(x){
#                   as.data.frame(z2[[x]]) %>%
#                     select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
#                            d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN) %>%
#                     mutate(ref = x, net_type = "pyr_only", model = "antag", int = "LLIN and IVM")
#                 }, simplify = F))
#

my_sim_antag_ITN_IVM <- function(){
  #pyr_out_list_antag_ITN_IVM <- purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
  pyr_out_list_antag_ITN_IVM <- lapply(pyr_param_list, antag_ITN_IVM_cov_loop)
  res_pyr_out_antag_ITN_IVM <- lapply(pyr_out_list_antag_ITN_IVM, runfun) #put these values into the model
  pyr_out_df_antag_ITN_IVM <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_ITN_IVM[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
  return(pyr_out_df_antag_ITN_IVM)

}

antag_ITN_IVM <- my_sim_antag_ITN_IVM()


#out <- out %>%
#  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
#                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
#                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
#                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
#                             TRUE ~ NA_character_))
#
#IVM_starts <- unique(out$IVRM_sr)
#MDA_6m_nets <- c(IVM_starts[2], IVM_starts[3], IVM_starts[4])
#MDA_1y_nets <- c(IVM_starts[5], IVM_starts[6], IVM_starts[7])
#MDA_2y_nets <- c(IVM_starts[8], IVM_starts[9], IVM_starts[10])
#
#out <- out %>%
#  mutate(net_age_MDA = "6 months")
#
#ggplot(out, aes(x = t, y= slide_prev0to5, col = as.factor(species)))+
#  geom_line()+
#  facet_wrap(vars(net_age_MDA))


##additive model
add_ITN_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  Q0_in <- itn_type_ivm_param[3]
  #IVRM_start_in <- itn_type_ivm_param[10] #failed because was getting made into numeric
  #IVRM_start_in <- ivm_nets_starting
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_constant_uptake.R", package = "ivRmectin"),
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
    #IVRM_start = IVRM_start_in,
    IVRM_start = ivm_parms1$IVRM_start,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_add_ITN <- function(){
  #pyr_out_list_add_ITN <- purrr::map2(y, x, add_ITN_cov_loop) #loop through all parameter values
  pyr_out_list_add_ITN <- lapply(pyr_param_list, add_ITN_cov_loop)
  res_pyr_out_add_ITN <- lapply(pyr_out_list_add_ITN, runfun) #put these values into the model
  pyr_out_df_add_ITN <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
    df <- as.data.frame(res_pyr_out_add_ITN[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
  return(pyr_out_df_add_ITN)

}

add_ITN <- my_sim_add_ITN()

#additive with IVM
add_ITN_IVM_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[7]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  Q0_in <- itn_type_ivm_param[3]
  #IVRM_start_in <- ivm_nets_starting
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_constant_uptake.R", package = "ivRmectin"),
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
    ivm_cov_par = 0.8, # proportion of popuulation receiving the endectocide
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

my_sim_add_ITN_IVM <- function(){
  #pyr_out_list_add_ITN_IVM <- purrr::map2(y, x, add_ITN_IVM_cov_loop) #loop through all parameter values
  pyr_out_list_add_ITN_IVM <- lapply(pyr_param_list, add_ITN_IVM_cov_loop)
  res_pyr_out_add_ITN_IVM <- lapply(pyr_out_list_add_ITN_IVM, runfun) #put these values into the model
  pyr_out_df_add_ITN_IVM <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
    df <- as.data.frame(res_pyr_out_add_ITN_IVM[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only"))}, simplify = F))
  return(pyr_out_df_add_ITN_IVM)

}

add_ITN_IVM <- my_sim_add_ITN_IVM()


#bound antag models
antag_LLIN <- antag_ITN %>%
  mutate(model = "antag_LLIN")
antag_LLIN_IVM <- antag_ITN_IVM %>%
  mutate(model = "antag_LLIN_IVM")
add_LLIN <- add_ITN %>%
  mutate(model = "add_LLIN")
add_LLIN_IVM <- add_ITN_IVM %>%
  mutate(model = "add_LLIN_IVM")

antag <- rbind(antag_LLIN, antag_LLIN_IVM)
add <- rbind(add_LLIN, add_LLIN_IVM)
saveRDS(antag, file = "analysis/exploring_interactions/MIM_poster/Q0/antag.rds")
saveRDS(add, file = "analysis/exploring_interactions/MIM_poster/Q0/add.rds")

#IVM only model and no intervention model (both can be model A)

#need to take out nets from pyr_param_df
pyr_param_df_no_nets <- pyr_param_df %>%
  select(init_EIR, Q0)
pyr_param_list_no_int <- list()
for (i in seq_len(nrow(pyr_param_df_no_nets))){
     pyr_param_list_no_int[[i]] <- as.numeric(pyr_param_df_no_nets[i,])
   }

#IVM only####
antag_IVM_cov_loop <- function(itn_type_ivm_param){
  init_EIR_in <- itn_type_ivm_param[1]
  Q0_in <- itn_type_ivm_param[2]
  #IVRM_start_in <- ivm_nets_starting
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    #num_int = 1,
    num_int = 1, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0.8, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_antag_IVM <- function(){
  #pyr_out_list_antag_ITN_IVM <- purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
  pyr_out_list_antag_IVM <- lapply(pyr_param_list_no_int, antag_IVM_cov_loop)
  res_pyr_out_antag_IVM <- lapply(pyr_out_list_antag_IVM, runfun) #put these values into the model
  pyr_out_df_antag_IVM <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
    df <- as.data.frame(res_pyr_out_antag_IVM[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only", model = "antag_IVM"))}, simplify = F))
  return(pyr_out_df_antag_IVM)

}

antag_IVM <- my_sim_antag_IVM()
saveRDS(antag_IVM, file = "analysis/exploring_interactions/MIM_poster/Q0/antag_IVM.rds")

#no intervention####
antag_no_int_cov_loop <- function(itn_type_ivm_param){
  init_EIR_in <- itn_type_ivm_param[1]
  Q0_in <- itn_type_ivm_param[2]
  #IVRM_start_in <- ivm_nets_starting
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    #num_int = 1,
    num_int = 1, # number of vector control (IRS and ITN) population groups
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
    Q0 = Q0_in
  )
  return(output)
}

my_sim_antag_no_int <- function(){
  #pyr_out_list_antag_ITN_IVM <- purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
  pyr_out_list_no_int<- lapply(pyr_param_list_no_int, antag_no_int_cov_loop)
  res_pyr_out_no_int <- lapply(pyr_out_list_no_int, runfun) #put these values into the model
  pyr_out_df_no_int <- do.call(rbind, sapply(1:(nrow(pyr_param_df)), function(x){
    df <- as.data.frame(res_pyr_out_no_int[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "pyr_only", model = "antag_no_int"))}, simplify = F))
  return(pyr_out_df_antag_IVM)

}

antag_no_int <- my_sim_antag_IVM()
saveRDS(antag_no_int, file = "analysis/exploring_interactions/MIM_poster/Q0/antag_no_int.rds")

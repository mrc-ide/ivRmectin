#extension of the ESA: looking at the impact of feeding behaviour and LLIN coverage on endpoints, with ivermectin MDA

require(tidyverse)
devtools::load_all()
require(tidyverse)

init_EIR <- 100

# Provide a value of the annual EIR for this model run
init_EIR <- 100 #low - 2, moderate - 15, high - 120

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years
time_period <- 365*10
IVM_begin <- (365*8)+180 #IVM MDA starts 2.5 y into the simulation
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int) #distributions take place every month for 3 months (3x300 micrograms/kg dosage)
ivm_cov = 0.8 #ivermectin coverage is set to 80% in these runs. assume that everyone is treated on the same day
#itn_cov_in = 0.8

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

itn_on <- 100 #nets introduced 100 days into the simulation --> reach equilibrium

net_seq <- seq(itn_on, 3650, by = 3*365) #net distributions are every 3 years

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

eff_len <- 23
ivm_on <- IVM_start[1]
ivm_off <- IVM_start[3]+eff_len

res_vector <- c(0, 0.1, 0.5, 0.7, 0.9)

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance == 0 | resistance  == 0.1 | resistance == 0.5 | resistance == 0.7 | resistance == 0.9)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

#IG2_d_ITN0 <- df_IG2$dn0_med
pyr_only_d_ITN0 <- df_pyr_only$dn0_med
#pyr_pbo_d_ITN0 <- df_pyr_pbo$dn0_med


pyr_param_list_bites_Bed <- list()
#pyr_pbo_param_list <- list()
#IG2_param_list <- list()


itn_cov_vector <- seq(0, 0.8, 0.2) #explore 0 to 80% resistance
#bites_Bed_vector <- seq(0.1, 0.9, 0.2)
#bites_Bed_vector <- seq(0.1, 0.9, 0.2)

#make a df for the species: here, we keep the Q0 of a gambiae-like vector and explore bites_Bed of many others
species <- c("gambiae", "arabiensis", "funestus", "stephensi")
#Q0 <- c("0.92", "0.8", "0.78", "0.51")
Q0 <- rep(0.92, 4)
#just going to change one thing at a time - keep Q0 high and explore multiple phi-Bs
bites_Bed <- c("0.85", "0.8", "0.78", "0.52")

biting_df <- data.frame(species, Q0, bites_Bed)

# for pyrethroid nets
pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vector, species = species)
pyr_param_df_1 <- left_join(pyr_param_df_crit,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df <- left_join(pyr_param_df_1,
                          biting_df %>% dplyr::select(species, Q0, bites_Bed),
                          by = c("species"))

#for other net brands
#pyr_pbo_param_df_crit <- expand.grid(dn0_med = pyr_pbo_d_ITN0, itn_cov = itn_cov_vector, species = species)
#pyr_pbo_param_df_1 <- left_join(pyr_pbo_param_df_crit,
#                            df_pyr_pbo %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
#                            by = c("dn0_med")) %>%
#  mutate(net_type = "pyrethroid PBO") %>%
#  mutate(gamman_med = gamman_med*365) %>%
#  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)
#
#pyr_pbo_param_df <- left_join(pyr_pbo_param_df_1,
#                          biting_df %>% dplyr::select(species, Q0, bites_Bed),
#                          by = c("species"))
#
#IG2_param_df_crit <- expand.grid(dn0_med = IG2_d_ITN0, itn_cov = itn_cov_vector, species = species)
#IG2_param_df_1 <- left_join(IG2_param_df_crit,
#                            df_IG2 %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
#                            by = c("dn0_med")) %>%
#  mutate(net_type = "IG2") %>%
#  mutate(gamman_med = gamman_med*365) %>%
#  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)
#
#IG2_param_df <- left_join(IG2_param_df_1,
#                          biting_df %>% dplyr::select(species, Q0, bites_Bed),
#                          by = c("species"))
#
#
#nets_param_df <- do.call("rbind", list(pyr_param_df, pyr_pbo_param_df,
#                                       IG2_param_df))
#
#
#pyr_param_df_arab <- pyr_param_df_gamb %>%
#  mutate(Q0 = 0.71,
#         bites_Bed = 0.8)

#pyr_param_df <- rbind(pyr_param_df_gamb, pyr_param_df_arab)

#no_int <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
#                                     num_int = 1,
#                                     #num_int = 2,
#                                     #ITN_IRS_on = 100,
#                                     #itn_cov = 0,
#                                     #het_brackets = 5,
#                                     #age = init_age,
#                                     init_EIR = 100,
#                                     #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
#                                     #admin2 = "Fatick",
#                                     ttt = ivm_parms1$ttt,
#                                     eff_len = ivm_parms1$eff_len,
#                                     haz = ivm_parms1$haz,
#                                     #ivm_cov_par = ivm_parms1$ivm_cov_par,
#                                     ivm_cov_par = 0,
#                                     ivm_min_age = ivm_parms1$ivm_min_age,
#                                     ivm_max_age = ivm_parms1$ivm_max_age,
#                                     IVRM_start = ivm_parms1$IVRM_start)
#
#no_int_run <- runfun(no_int)
#df_no_int_run <- as.data.frame(no_int_run)

for (i in seq_len(nrow(pyr_param_df))){
  pyr_param_list_bites_Bed[[i]] <- as.numeric(pyr_param_df[i,])
}

#for (i in seq_len(nrow(pyr_pbo_param_df))){
#  pyr_pbo_param_list[[i]] <- as.numeric(pyr_pbo_param_df[i,])
#}
#
#for (i in seq_len(nrow(IG2_param_df))){
#  IG2_param_list[[i]] <- as.numeric(IG2_param_df[i,])
#}


#
#no_res_d_ITN0 <- pyr_param_df$d_ITN0[1]
#no_res_r_ITN0 <- pyr_param_df$r_ITN0[1]
#no_res_itn_half_life <- pyr_param_df$itn_half_life[1]
#

antag_ITN_IVM_loop_bites_Bed <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in<- itn_type_ivm_param[9]
  Q0_in <- itn_type_ivm_param[8]
  r_ITN0_in <- itn_type_ivm_param[4]
  itn_half_life_in <- itn_type_ivm_param[5]
  #itn_half_life_in <- itn_type_ivm_param[4]*365
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = itn_on,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0.8, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in,
    itn_half_life = itn_half_life_in
  )
  return(output)
}


pyr_out_list_bites_Bed <- lapply(pyr_param_list_bites_Bed, antag_ITN_IVM_loop_bites_Bed) #loop through all parameter values

res_pyr_out_bites_Bed <- lapply(pyr_out_list_bites_Bed, runfun) #put these values into the model

bites_Bed_vec <- bites_Bed

antag_pyr_out_df_bites_Bed <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)), function(x){
                                  as.data.frame(res_pyr_out_bites_Bed[[x]]) %>%
                                    select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, bites_Bed,
                                           d_ITN0, r_ITN0, itn_loss, Q0) %>%
                                    mutate(ref = x, net_type = "pyr_only",
                                           species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
                                                               bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
                                                               bites_Bed == bites_Bed_vec[3] ~ "funestus",
                                                               bites_Bed == bites_Bed_vec[4] ~ "stephensi",
                                                               TRUE ~ NA_character_))
                                }, simplify = F))

#running separately for each net type, to help with speed and checks

#pyr_pbo_out_list <- lapply(pyr_pbo_param_list, antag_ITN_cov_loop) #loop through all parameter values
#
#res_pyr_pbo_out <- lapply(pyr_pbo_out_list, runfun) #put these values into the model
#
#pyr_pbo_out_df <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)), function(x){
#  as.data.frame(res_pyr_out[[x]]) %>%
#    select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, bites_Bed,
#           d_ITN0, r_ITN0, itn_loss, Q0) %>%
#    mutate(ref = x, net_type = "pyr_pbo",
#           species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
#                               bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
#                               bites_Bed == bites_Bed_vec[3] ~ "funestus",
#                               bites_Bed == bites_Bed_vec[4] ~ "stephensi",
#                               TRUE ~ NA_character_))
#}, simplify = F))
#
#tail(pyr_pbo_out_df)
#
#
#IG2_out_list <- lapply(IG2_param_list, antag_ITN_cov_loop) #loop through all parameter values
#
#res_IG2_out <- lapply(IG2_out_list, runfun) #put these values into the model
#
#IG2_out_df <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)), function(x){
#  as.data.frame(res_pyr_out[[x]]) %>%
#    select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, bites_Bed,
#           d_ITN0, r_ITN0, itn_loss, Q0) %>%
#    mutate(ref = x, net_type = "IG2",
#           species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
#                               bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
#                               bites_Bed == bites_Bed_vec[3] ~ "funestus",
#                               bites_Bed == bites_Bed_vec[4] ~ "stephensi",
#                               TRUE ~ NA_character_))
#}, simplify = F))
#
#head(IG2_out_df)
#
##write out files

write.csv(antag_pyr_out_df_bites_Bed, file = "analysis/biting_LLINs/bites_Bed/biting_antag_out_pyr_IVM_ITN_bites_Bed_df.csv", row.names = FALSE)
#write.csv(pyr_pbo_out_df, file = "analysis/biting_LLINs/biting_antag_out_pyr_pbo_df.csv", row.names = FALSE)
#write.csv(IG2_out_df, file = "analysis/biting_LLINs/biting_antag_out_IG2_df.csv", row.names = FALSE)

#additive model
add_ITN_IVM_loop_bites_Bed <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in<- itn_type_ivm_param[9]
  Q0_in <- itn_type_ivm_param[8]
  r_ITN0_in <- itn_type_ivm_param[4]
  itn_half_life_in <- itn_type_ivm_param[5]
  #itn_half_life_in <- itn_type_ivm_param[4]*365
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide_constant_uptake.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = itn_on,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0.8, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in,
    itn_half_life = itn_half_life_in
  )
  return(output)
}

add_pyr_out_list_bites_Bed <- lapply(pyr_param_list_bites_Bed, add_ITN_IVM_loop_bites_Bed) #loop through all parameter values

add_res_pyr_out_bites_Bed <- lapply(add_pyr_out_list_bites_Bed, runfun) #put these values into the model


add_pyr_out_df_bites_Bed <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)), function(x){
  as.data.frame(add_res_pyr_out_bites_Bed[[x]]) %>%
    select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, bites_Bed,
           d_ITN0, r_ITN0, itn_loss, Q0) %>%
    mutate(ref = x, net_type = "pyr_only",
           species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
                               bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
                               bites_Bed == bites_Bed_vec[3] ~ "funestus",
                               bites_Bed == bites_Bed_vec[4] ~ "stephensi",
                               TRUE ~ NA_character_))
}, simplify = F))

write.csv(add_pyr_out_df_bites_Bed, file = "analysis/biting_LLINs/bites_Bed/biting_add_out_pyr_IVM_ITN_bites_Bed_df.csv", row.names = FALSE)


#running separately for each net type, to help with speed and checks

#add_pyr_pbo_out_list <- lapply(pyr_pbo_param_list, add_ITN_cov_loop) #loop through all parameter values
#
#add_res_pyr_pbo_out <- lapply(add_pyr_pbo_out_list, runfun) #put these values into the model
#
#add_pyr_pbo_out_df <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)), function(x){
#  as.data.frame(add_res_pyr_pbo_out[[x]]) %>%
#    select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, bites_Bed,
#           d_ITN0, r_ITN0, itn_loss, Q0) %>%
#    mutate(ref = x, net_type = "pyr_pbo",
#           species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
#                               bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
#                               bites_Bed == bites_Bed_vec[3] ~ "funestus",
#                               bites_Bed == bites_Bed_vec[4] ~ "stephensi",
#                               TRUE ~ NA_character_))
#}, simplify = F))
#
#tail(add_pyr_pbo_out_df)
#
#
#add_IG2_out_list <- lapply(IG2_param_list, antag_ITN_cov_loop) #loop through all parameter values
#
#add_res_IG2_out <- lapply(add_IG2_out_list, runfun) #put these values into the model
#
#add_IG2_out_df <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)), function(x){
#  as.data.frame(add_res_IG2_out[[x]]) %>%
#    select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, bites_Bed,
#           d_ITN0, r_ITN0, itn_loss, Q0) %>%
#    mutate(ref = x, net_type = "IG2",
#           species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
#                               bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
#                               bites_Bed == bites_Bed_vec[3] ~ "funestus",
#                               bites_Bed == bites_Bed_vec[4] ~ "stephensi",
#                               TRUE ~ NA_character_))
#}, simplify = F))
#
#head(add_IG2_out_df)


#BASELINE MODELS: LLIN ONLY####

#baseline - LLIN only, antagonistic, phi-B#
antag_ITN_cov_loop_bites_Bed <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in<- itn_type_ivm_param[9]
  Q0_in <- itn_type_ivm_param[8]
  r_ITN0_in <- itn_type_ivm_param[4]
  itn_half_life_in <- itn_type_ivm_param[5]
  #itn_half_life_in <- itn_type_ivm_param[4]*365
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = itn_on,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide. NO IVERMECTIN
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

antag_pyr_out_list_LLIN_bites_Bed <- lapply(pyr_param_list_bites_Bed, antag_ITN_cov_loop_bites_Bed) #loop through all parameter values

res_antag_pyr_out_LLIN_bites_Bed <- lapply(antag_pyr_out_list_LLIN_bites_Bed, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_antag_LLIN_bites_Bed <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)),
                                               function(x){
                              as.data.frame(res_antag_pyr_out_LLIN_bites_Bed[[x]]) %>%
                                select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0) %>%
                                mutate(ref = x, net_type = "pyr_only", model = "antag",
                                       species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
                                                           bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
                                                           bites_Bed == bites_Bed_vec[3] ~ "funestus",
                                                           bites_Bed == bites_Bed_vec[4] ~ "stephensi",
                                                           TRUE ~ NA_character_))
                            }, simplify = F))

write.csv(pyr_out_df_antag_LLIN_bites_Bed, file = "analysis/biting_LLINs/bites_Bed/pyr_out_df_antag_LLIN_bites_Bed.csv", row.names = FALSE)


#baseline model, additive

#baseline - LLIN only, antagonistic, phi-B#
add_ITN_cov_loop_bites_Bed <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in<- itn_type_ivm_param[9]
  Q0_in <- itn_type_ivm_param[8]
  r_ITN0_in <- itn_type_ivm_param[4]
  itn_half_life_in <- itn_type_ivm_param[5]
  #itn_half_life_in <- itn_type_ivm_param[4]*365
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide_constant_uptake.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = itn_on,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

add_pyr_out_list_LLIN_bites_Bed <- lapply(pyr_param_list_bites_Bed, add_ITN_cov_loop_bites_Bed) #loop through all parameter values

res_add_pyr_out_LLIN_bites_Bed <- lapply(add_pyr_out_list_LLIN_bites_Bed, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_add_LLIN_bites_Bed <- do.call(rbind, sapply(1:(length(itn_cov_vector)*length(species)*length(res_vector)),
                                                         function(x){
                                                           as.data.frame(res_add_pyr_out_LLIN_bites_Bed[[x]]) %>%
                                                             select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                                                    d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0) %>%
                                                             mutate(ref = x, net_type = "pyr_only", model = "antag",
                                                                    species = case_when(bites_Bed == bites_Bed_vec[1] ~ "gambiae",
                                                                                        bites_Bed == bites_Bed_vec[2] ~ "arabiensis",
                                                                                        bites_Bed == bites_Bed_vec[3] ~ "funestus",
                                                                                        bites_Bed == bites_Bed_vec[4] ~ "stephensi",
                                                                                        TRUE ~ NA_character_))
                                                         }, simplify = F))

write.csv(pyr_out_df_add_LLIN_bites_Bed, file = "analysis/biting_LLINs/bites_Bed/pyr_out_df_add_LLIN_bites_Bed.csv", row.names = FALSE)

#model runs for the exploring interactions paper

#going to look at dynamics in EIR, prevalence for vectors in different transmission settings

#the vector parameters are for: gambiae, arabiensis, funestus, stephensi
#first 4 values in vector are to explore parameter space


#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
require(tidyverse)
bites_Bed_vector <- c(0.85, 0.8, 0.78, 0.52)
Q0_vector <-       c(0.92, 0.71, 0.94, 0.21)

species_df <- data.frame(bites_Bed_vector, Q0_vector)
species_df$species <- c("gambiae", "arabiensis", "funestus", "stephensi")

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

#ivm on when nets are 6 months old
IVM_begin1 <- (365*6)+180 # 6 months into new net distribution
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

#when nets are 1yo

IVM_begin2 <- 365*7
IVM_start2 <- c(IVM_begin2, IVM_begin2+mda_int, IVM_begin2 + mda_int + mda_int)


#when nets are 2.5yo
IVM_begin3 <- (365*8)+180
IVM_start3 <- c(IVM_begin3, IVM_begin3+mda_int, IVM_begin3+mda_int+mda_int)


ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

itn_on <- 100 #introduce nets 100 days into simulation

net_seq <- seq(100, 3650, by = 3*365)

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#IVM_starters <- c(IVM_start1[1], IVM_start2[1], IVM_start3[1])

IVM_start <- numeric(length = 3)

IVM_starting <- list(IVM_start1, IVM_start2, IVM_start3)

#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivm_fun(
  #IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)


ivm_parms2 <- function(IVM_starting) {
  ivm_fun(
    #IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
    IVM_start_times = IVM_starting,
    time_period = time_period,
    hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
    #hazard_profile = hazzy,
    ivm_coverage=ivm_cov,
    ivm_min_age=5,
    ivm_max_age = 90)
}

ivm_parmscheck <- ivm_fun(
  IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  #IVM_start_times = c(2370, 2400, 2430),
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms3 <- lapply(IVM_starting, ivm_parms2)

#extract the ivermectin info from this
ivm_nets_1 <- ivm_parms3[[1]]$IVRM_start #6m nets
ivm_nets_2 <- ivm_parms3[[2]]$IVRM_start #1y nets
ivm_nets_3 <- ivm_parms3[[3]]$IVRM_start #2.5y nets

ivm_nets_starting <- list(ivm_nets_1, ivm_nets_2, ivm_nets_3)


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
pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vector, bites_Bed = bites_Bed_vector,
                                 init_EIR = init_EIR_vec)
pyr_param_df <- left_join(pyr_param_df_crit,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med) %>%
  left_join(species_df, by = c("bites_Bed" = "bites_Bed_vector")) %>%
  rename(Q0 = Q0_vector)


#pyr_param_df2 <- crossing(pyr_param_df, ivm_nets_starting) %>%
#  mutate(scenario_number = row_number())
saveRDS(pyr_param_df, file = "data/species_runs_net_age.rds")

#saveRDS(pyr_param_df2, file = "data/species_runs_net_age.rds")


for (i in seq_len(nrow(pyr_param_df))){
  pyr_param_list[[i]] <- as.numeric(pyr_param_df[i,])
}

#for (i in seq_len(nrow(pyr_param_df2))){
#  pyr_param_list[[i]] <- as.list(pyr_param_df2[i,])
#}

#ivm_input <- numeric(length = 3651)
#going to pass in the ivm_nets_starting list and access things
antag_ITN_cov_loop <- function(itn_type_ivm_param, ivm_nets_starting){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[8]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  #IVRM_start_in <- itn_type_ivm_param[10] #failed because was getting made into numeric
  IVRM_start_in <- ivm_nets_starting
  #itn_half_life_in <- itn_type_ivm_param[4]*365
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
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
    IVRM_start = IVRM_start_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

#pyr_out_list_antag_ITN <- lapply(pyr_param_list, ivm_nets_starting = ivm_nets_starting, antag_ITN_cov_loop) #loop through all parameter values
x <- rep(ivm_nets_starting, nrow(pyr_param_df)) #do reps to ensure have the same length
#need to rep the list by 3
pyr_out_list_antag_ITN <- purrr::map2(pyr_param_list, x, antag_ITN_cov_loop) #loop through all parameter values
res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_antag_ITN <- do.call(rbind,
                                sapply(1:(nrow(pyr_param_df*3)), function(x){
                                              as.data.frame(res_pyr_out_antag_ITN[[x]]) %>%
                                                select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN) %>%
                                                mutate(ref = x, net_type = "pyr_only", model = "antag", int = "LLIN")
                                            }, simplify = F))

#write.csv(pyr_out_df_antag_ITN, file = "analysis/exploring_interactions/model_output/antag_pyr_LLIN_nets_age.csv", row.names = FALSE)
saveRDS(pyr_out_df_antag_ITN, file = "analysis/exploring_interactions/model_output/antag_pyr_LLIN_nets_age.rds")
#antag model with ivermectin and nets

antag_ITN_IVM_cov_loop <- function(itn_type_ivm_param, ivm_nets_starting){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[8]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  IVRM_start_in <- ivm_nets_starting
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
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
    IVRM_start = IVRM_start_in ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

#pyr_out_list_antag_ITN_IVM <- lapply(pyr_param_list, ivm_nets_starting = ivm_nets_starting[1], antag_ITN_IVM_cov_loop) #loop through all parameter values
###testing###
x <- rep(ivm_nets_starting, 2) #do reps to ensure have the same length
y <- rep(pyr_param_list[1:2], 3)
#
z <-purrr::map2(y, x, antag_ITN_IVM_cov_loop) #loop through all parameter values
z2 <- lapply(z, runfun) #put these values into the model
z2df  <- do.call(rbind,
                                         sapply(1:(6), function(x){
                                           as.data.frame(z2[[x]]) %>%
                                             select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                                    d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN) %>%
                                             mutate(ref = x, net_type = "pyr_only", model = "antag", int = "LLIN and IVM")
                                         }, simplify = F))

ggplot(z2df, aes(x = t, y = EIR_tot))+
    geom_line()



pyr_out_list_antag_ITN_IVM <-purrr::map2(pyr_param_list, x, antag_ITN_IVM_cov_loop) #loop through all parameter values

res_pyr_out_antag_ITN_IVM <- lapply(pyr_out_list_antag_ITN_IVM, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_antag_ITN_IVM <- do.call(rbind,
                                    sapply(1:(nrow(pyr_param_df*3)), function(x){
                                                as.data.frame(res_pyr_out_antag_ITN_IVM[[x]]) %>%
                                                  select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                                         d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN) %>%
                                                  mutate(ref = x, net_type = "pyr_only", model = "antag", int = "LLIN and IVM")
                                              }, simplify = F))
saveRDS(pyr_out_df_antag_ITN_IVM, file = "analysis/exploring_interactions/model_output/antag_pyr_LLIN_IVM_nets_age.rds")

#write.csv(pyr_out_df_antag_ITN_IVM, file = "analysis/exploring_interactions/model_output/antag_pyr_LLIN_IVM_nets_age.csv", row.names = FALSE)

#additive, nets only
add_ITN_cov_loop <- function(itn_type_ivm_param, ivm_nets_starting){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[8]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  IVRM_start_in <- ivm_nets_starting
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide_constant_uptake.R",
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
    IVRM_start = IVRM_start_in ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

pyr_out_list_add_ITN <- purrr::map2(pyr_param_list, x,add_ITN_cov_loop) #loop through all parameter values

res_pyr_out_add_ITN <- lapply(pyr_out_list_add_ITN, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_add_ITN <- do.call(rbind,
                              sapply(1:(nrow(pyr_param_df*3)), function(x){
                                          as.data.frame(res_pyr_out_add_ITN[[x]]) %>%
                                            select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                                   d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN) %>%
                                            mutate(ref = x, net_type = "pyr_only", model = "add", int = "LLIN")
                                        }, simplify = F))

saveRDS(pyr_out_df_add_ITN, file = "analysis/exploring_interactions/model_output/add_pyr_LLIN.csv")
#write.csv(pyr_out_df_add_ITN, file = "analysis/exploring_interactions/model_output/add_pyr_LLIN.csv", row.names = FALSE)

#add: LLIN and IVM
add_ITN_IVM_cov_loop <- function(itn_type_ivm_param, ivm_nets_starting){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[8]
  init_EIR_in <- itn_type_ivm_param[4]
  r_ITN0_in <- itn_type_ivm_param[5]
  itn_half_life_in <- itn_type_ivm_param[6]
  IVRM_start_in <- ivm_nets_starting
  #itn_half_life_in <- itn_type_ivm_param[4]*365
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide_constant_uptake.R",
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
    IVRM_start = IVRM_start_in ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

pyr_out_list_add_ITN_IVM <- purrr::map2(pyr_param_list, x, add_ITN_IVM_cov_loop) #loop through all parameter values

res_pyr_out_add_ITN_IVM <- lapply(pyr_out_list_add_ITN_IVM, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_add_ITN_IVM <- do.call(rbind,
                                  sapply(1:(nrow(pyr_param_df*3)), function(x){
                                              as.data.frame(res_pyr_out_add_ITN_IVM[[x]]) %>%
                                                select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0, IVRM_sr, s_ITN, d_ITN, r_ITN) %>%
                                                mutate(ref = x, net_type = "pyr_only", model = "add", int = "LLIN and IVM")
                                            }, simplify = F))
saveRDS(pyr_out_df_add_ITN_IVM, file = "analysis/exploring_interactions/model_output/add_pyr_LLIN_IVM.rds")
#write.csv(pyr_out_df_add_ITN_IVM, file = "analysis/exploring_interactions/model_output/add_pyr_LLIN_IVM.csv", row.names = FALSE)


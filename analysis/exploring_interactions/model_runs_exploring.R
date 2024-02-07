#model runs for the exploring interactions paper

#going to look at dynamics in EIR, prevalence for vectors in different transmission settings

#the vector parameters are for: gambiae, arabiensis, funestus, stephensi
#first 4 values in vector are to explore parameter space
bites_Bed_vector <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.85, 0.8, 0.78, 0.52)
Q0_vector <-       c(0.1, 0.3, 0.5, 0.7, 0.9, 0.92, 0.71, 0.94, 0.21)

#4 figures:
#1) Dynamics plot for A.gambiae-like vector at 10% resistance in different transmission settings
#2) Efficacy plot: predicted by anatagonisitic and additive model, shapes and colours for different phi-B, facet by Q0 and transmission setting
#3) Efficacy by species: creat a geom tile, for each transmission setting, and show the relative different in EIR and prevalence. Y is resistance, x is species

#relatives are LLIN & IVM compared to LLIN only
#ivRmectin model has 0.89 as bites_Bed gambiae, going to reset as 0.85 from PNAS paper

require(tidyverse)

#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
require(tidyverse)
# Create a vector of age categories for the model
#init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR_vec <- c(2, 25, 100) #low - 2, moderate - 15, high - 120 --> Ellie: low = 2, med = 25, high = 100

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years
time_period <- 365*10
IVM_begin <- (365*8)+180
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)
ivm_cov = 0.8
itn_cov_in = 0.8

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

pyr_param_list <- list()

#repeat for pyrethroid nets
pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vector, bites_Bed = bites_Bed_vector,
                                 Q0 = Q0_vector, init_EIR = init_EIR_vec)
pyr_param_df <- left_join(pyr_param_df_crit,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

for (i in seq_len(nrow(pyr_param_df))){
  pyr_param_list[[i]] <- as.numeric(pyr_param_df[i,])
}

antag_ITN_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[4]
  init_EIR_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
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
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

pyr_out_list_antag_ITN <- lapply(pyr_param_list, antag_ITN_cov_loop) #loop through all parameter values

res_pyr_out_antag_ITN <- lapply(pyr_out_list_antag_ITN, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_antag_ITN <- do.call(rbind,
                            sapply(1:(length(itn_cov_vector)*length(res_vector)* length(bites_Bed_vector)*
                                        length(Q0_vector)*length(init_EIR_vec)), function(x){
                              as.data.frame(res_pyr_out_antag_ITN[[x]]) %>%
                                select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                       d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0) %>%
                                mutate(ref = x, net_type = "pyr_only", model = "antag", int = "LLIN")
                            }, simplify = F))

write.csv(pyr_out_df_antag_ITN, file = "analysis/exploring_interactions/model_output/antag_pyr_LLIN.csv", row.names = FALSE)

#antag model with ivermectin and nets

antag_ITN_IVM_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[4]
  init_EIR_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
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
    ivm_cov_par = 0.8, # proportion of popuulation receiving the endectocide
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

pyr_out_list_antag_ITN_IVM <- lapply(pyr_param_list, antag_ITN_IVM_cov_loop) #loop through all parameter values

res_pyr_out_antag_ITN_IVM <- lapply(pyr_out_list_antag_ITN_IVM, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_antag_ITN_IVM <- do.call(rbind,
                                sapply(1:(length(itn_cov_vector)*length(res_vector)* length(bites_Bed_vector)* length(Q0_vector)
                                          *length(init_EIR_vec)), function(x){
                                  as.data.frame(res_pyr_out_antag_ITN_IVM[[x]]) %>%
                                    select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                           d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0) %>%
                                    mutate(ref = x, net_type = "pyr_only", model = "antag", int = "LLIN and IVM")
                                }, simplify = F))

write.csv(pyr_out_df_antag_ITN_IVM, file = "analysis/exploring_interactions/model_output/antag_pyr_LLIN_IVM.csv", row.names = FALSE)

#additive, nets only
add_ITN_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[4]
  init_EIR_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
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

pyr_out_list_add_ITN <- lapply(pyr_param_list, add_ITN_cov_loop) #loop through all parameter values

res_pyr_out_add_ITN <- lapply(pyr_out_list_add_ITN, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_add_ITN <- do.call(rbind,
                                    sapply(1:(length(itn_cov_vector)*length(res_vector)* length(bites_Bed_vector)* length(Q0_vector)
                                              *length(init_EIR_vec)), function(x){
                                      as.data.frame(res_pyr_out_add_ITN[[x]]) %>%
                                        select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                               d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0) %>%
                                        mutate(ref = x, net_type = "pyr_only", model = "add", int = "LLIN")
                                    }, simplify = F))

write.csv(pyr_out_df_add_ITN, file = "analysis/exploring_interactions/model_output/add_pyr_LLIN.csv", row.names = FALSE)

#add: LLIN and IVM
add_ITN_IVM_cov_loop <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  Q0_in <- itn_type_ivm_param[4]
  init_EIR_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
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
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    bites_Bed = bites_Bed_in,
    Q0 = Q0_in
  )
  return(output)
}

pyr_out_list_add_ITN_IVM <- lapply(pyr_param_list, add_ITN_IVM_cov_loop) #loop through all parameter values

res_pyr_out_add_ITN_IVM <- lapply(pyr_out_list_add_ITN_IVM, runfun) #put these values into the model

#bites_Bed_vec <- c(0.85, 0.8)

pyr_out_df_add_ITN_IVM <- do.call(rbind,
                              sapply(1:(length(itn_cov_vector)*length(res_vector)* length(bites_Bed_vector)* length(Q0_vector)
                                        *length(init_EIR_vec)), function(x){
                                          as.data.frame(res_pyr_out_add_ITN_IVM[[x]]) %>%
                                            select(t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5,
                                                   d_ITN0, r_ITN0, itn_loss, bites_Bed, Q0) %>%
                                            mutate(ref = x, net_type = "pyr_only", model = "add", int = "LLIN and IVM")
                                        }, simplify = F))

write.csv(pyr_out_df_add_ITN, file = "analysis/exploring_interactions/model_output/add_pyr_LLIN_IVM.csv", row.names = FALSE)


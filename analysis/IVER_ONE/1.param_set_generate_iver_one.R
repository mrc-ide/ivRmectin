#getting parameter set for iver-one

#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)

#KWALE Q0 from Marta

#funestus: 76 of 107 bloodmeals on humans
#arabiensis: 3 of 13 bloodmeals on humans
#rivolurum: 6 of 9 bloodmeals on humans

prop_funestus <- 106/(106+13+9+12)
prop_arab <- 13/(106+13+9+12)
prop_rivo <- 9/(106+13+9+12)
prop_gamb <- 12/(106+13+9+12)

hbi_mosq <- (((76/106)*prop_funestus) + ((3/13)*prop_arab) + ((3/9)*prop_rivo)) + ((5/12)*prop_gamb)
hbi_mosq <- round(hbi_mosq, 2)
Q0_kwale <- hbi_mosq


#weighted estimate from Ellie's site file
bites_Bed_kwale <- (0.85*0.33) + (0.8*0.45) + (0.78*0.22)
itn_cov_kwale <- 0.58 #average from site file since 2018
res_kwale <- 0.27

init_EIR_kwale <- 30

ivm_cov <- 0.7

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance %in% res_kwale)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

d_ITN0_kwale <- df_pyr_only$dn0_med


#pyr_param_list <- list()

pyr_param_df_crit_kwale <- expand.grid(dn0_med = d_ITN0_kwale, itn_cov = itn_cov_kwale,
                                 init_EIR = init_EIR_kwale,
                                 bites_Bed = bites_Bed_kwale, ivm_cov = ivm_cov,
                                 Q0 = Q0_kwale)
pyr_param_df_kwale <- left_join(pyr_param_df_crit_kwale,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df_kwale <- pyr_param_df_kwale %>%
  select(-resistance) #remove the resistance column.

saveRDS(pyr_param_df_kwale, "analysis/IVER_ONE/kwale_param.rds")

#BUSIA
Q0_busia <- 0.83
bites_Bed_busia <- 0.82
itn_cov_busia <- 0.63
res_busia <- 0.44
init_EIR_busia <- 30 #update later

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance %in% res_busia)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

d_ITN0_busia <- df_pyr_only$dn0_med


pyr_param_df_crit_busia <- expand.grid(dn0_med = d_ITN0_busia, itn_cov = itn_cov_busia,
                                       init_EIR = init_EIR_busia,
                                       bites_Bed = bites_Bed_busia, ivm_cov = ivm_cov,
                                       Q0 = Q0_busia)

pyr_param_df_busia <- left_join(pyr_param_df_crit_busia,
                                df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                                by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df_busia <- pyr_param_df_busia %>%
  select(-resistance) #remove the resistance column.

saveRDS(pyr_param_df_busia, "analysis/IVER_ONE/busia_param.rds")

#HOMA BAY
Q0_homa_bay <- 0.81
bites_Bed_homa_bay <- 0.81
itn_cov_homa_bay <- 0.57
res_homa_bay <- 0.36
init_EIR_homa_bay <- 30 #update later. LOW SETTING

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance %in% res_homa_bay)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

d_ITN0_homa_bay <- df_pyr_only$dn0_med


pyr_param_df_crit_homa_bay <- expand.grid(dn0_med = d_ITN0_homa_bay, itn_cov = itn_cov_homa_bay,
                                       init_EIR = init_EIR_homa_bay,
                                       bites_Bed = bites_Bed_homa_bay, ivm_cov = ivm_cov,
                                       Q0 = Q0_homa_bay)

pyr_param_df_homa_bay <- left_join(pyr_param_df_crit_homa_bay,
                                df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                                by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df_homa_bay <- pyr_param_df_homa_bay %>%
  select(-resistance) #remove the resistance column.

saveRDS(pyr_param_df_homa_bay, "analysis/IVER_ONE/homa_bay_param.rds")

#MIGORI
Q0_migori <- 0.82
bites_Bed_migori <- 0.82
itn_cov_migori <- 0.57
res_migori <- 0.38
init_EIR_migori <- 30 #update later. LOW SETTING

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance %in% res_migori)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

d_ITN0_migori <- df_pyr_only$dn0_med


pyr_param_df_crit_migori <- expand.grid(dn0_med = d_ITN0_migori, itn_cov = itn_cov_migori,
                                          init_EIR = init_EIR_migori,
                                          bites_Bed = bites_Bed_migori, ivm_cov = ivm_cov,
                                          Q0 = Q0_migori)

pyr_param_df_migori <- left_join(pyr_param_df_crit_migori,
                                   df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                                   by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df_migori <- pyr_param_df_migori %>%
  select(-resistance) #remove the resistance column.

saveRDS(pyr_param_df_migori, "analysis/IVER_ONE/migori_param.rds")

#SIAYA
Q0_siaya <- 0.82
bites_Bed_siaya <- 0.82
itn_cov_siaya <- 0.61
res_siaya <- 0.49
init_EIR_siaya <- 30 #update later. LOW SETTING

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance %in% res_siaya)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

d_ITN0_siaya <- df_pyr_only$dn0_med


pyr_param_df_crit_siaya <- expand.grid(dn0_med = d_ITN0_siaya, itn_cov = itn_cov_siaya,
                                        init_EIR = init_EIR_siaya,
                                        bites_Bed = bites_Bed_siaya, ivm_cov = ivm_cov,
                                        Q0 = Q0_siaya)

pyr_param_df_siaya <- left_join(pyr_param_df_crit_siaya,
                                 df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                                 by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df_siaya <- pyr_param_df_siaya %>%
  select(-resistance) #remove the resistance column.

saveRDS(pyr_param_df_siaya, "analysis/IVER_ONE/siaya_param.rds")






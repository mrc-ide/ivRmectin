#getting parameter set for Bijagos scenario (moderate pyrethroid resistance) that will be used for malariasim validation

#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
#bites_Bed_vector_sens <- c(0.25, 0.5, 0.75, 0.9)
#Q0_vector_sens <- c(0.25, 0.5, 0.75, 0.9)
#llin_cov_vector_sens <- c(0, 0.2, 0.6, 0.8)
#res_vector_sens <- c(0, 0.2, 0.6, 0.8)
phi_vals <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/field work/fieldwork-data/cleaned/phi_calc/phi_calc_stan.csv")
phi_vals
bites_Bed_vec <-c(min(phi_vals$phi_lower_round), max(phi_vals$phi_upper_round)) #from my Bijagos paper
bites_Bed_vec <- c(0.1, 0.9) #override to explore extremes
Q0_vec <- 0.92
Q0_vec <- c(0.1, 0.9) #override to explore extremes
res <- c(0.5, 0.9)
init_EIR_vec <- c(30, 100)
itn_cov_vec <- c(0.5, 0.9) # two extremes
ivm_cov_vec <- c(0.5, 0.9) # two extremes

path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
filenames <- list.files(path = path_nets, pattern = ".csv$", full.names = TRUE)
lapply(filenames, function(x) {
  df <- subset(read.csv(x), resistance %in% res)
  return(df)
}) -> list_data

names(list_data) <- c("df_pyr_only", "df_pyr_pbo", "df_IG2")
list2env(list_data, .GlobalEnv)

pyr_only_d_ITN0 <- df_pyr_only$dn0_med

pyr_param_list_og <- list()
pyr_param_list <- list()

#repeat for pyrethroid nets
pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vec,
                                 init_EIR = init_EIR_vec,
                                 bites_Bed = bites_Bed_vec, ivm_cov = ivm_cov_vec,
                                 Q0 = Q0_vec)
pyr_param_df <- left_join(pyr_param_df_crit,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df <- pyr_param_df %>%
  select(-resistance) #remove the resistance column. We know corresponds to 0.5 and 0.9

saveRDS(pyr_param_df, "analysis/exploring_interactions/malariasim-odin/scenario-parameter-set.rds")
#and save it to the DIDE drive.
saveRDS(pyr_param_df, "W:/endectocides-cluster/data/scenario-parameter-set.rds")

pyr_param_df_no_int <- pyr_param_df %>%
  select(init_EIR, ivm_cov, Q0)

saveRDS(pyr_param_df_no_int, "analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-no-int.rds")
#and save it to the DIDE drive.
saveRDS(pyr_param_df_no_int, "W:/endectocides-cluster/data/scenario-parameter-set-no-int.rds")


#same but with different endec_mu and wane
endec_mu_vec <- seq(0, 1, 0.001)
wane_vec <- seq(0,0.1,0.001)

pyr_param_df_crit2 <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov_vec,
                                 init_EIR = init_EIR_vec,
                                 bites_Bed = bites_Bed_vec, Q0 = Q0_vec, ivm_cov = ivm_cov_vec, endec_mu = endec_mu_vec,
                                 wane = wane_vec)
pyr_param_df2 <- left_join(pyr_param_df_crit2,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

saveRDS(pyr_param_df2, "analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-exp-decay.rds")
saveRDS(pyr_param_df2, "W:/endectocides-cluster/data/scenario-parameter-set-exp-decay.rds")

pyr_param_df2_no_int <- pyr_param_df2 %>%
  select(init_EIR, Q0, endec_mu, wane)

saveRDS(pyr_param_df2_no_int, "analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-exp-decay-no-int.rds")
saveRDS(pyr_param_df2_no_int, "W:/endectocides-cluster/data/scenario-parameter-set-exp-decay-no-int.rds")

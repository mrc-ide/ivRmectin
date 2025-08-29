#getting parameter set for simple MDA impact in relation to timing of ITN campaign.

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
bites_Bed <- 0.95
#Q0_vec <- 0.92
Q0 <- 0.95
res <- 0.55
init_EIR_vec <- c(30, 100)
itn_cov <- 0.5
ivm_cov <- 0.7

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
pyr_param_df_crit <- expand.grid(dn0_med = pyr_only_d_ITN0, itn_cov = itn_cov,
                                 init_EIR = init_EIR_vec,
                                 bites_Bed = bites_Bed, ivm_cov = ivm_cov,
                                 Q0 = Q0)
pyr_param_df <- left_join(pyr_param_df_crit,
                          df_pyr_only %>% dplyr::select(dn0_med, rn0_med, gamman_med, resistance),
                          by = c("dn0_med")) %>%
  #mutate(net_type = "pyrethroid only") %>%
  mutate(gamman_med = gamman_med*365) %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med, itn_half_life = gamman_med)

pyr_param_df <- pyr_param_df %>%
  select(-resistance) #remove the resistance column. We know corresponds to 0.5 and 0.9

saveRDS(pyr_param_df, "analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-mda-times.rds")






























#more tricky when nets on because mu won't be 0.132


#and save it to the DIDE drive.
#saveRDS(pyr_param_df, "W:/endectocides-cluster/data/scenario-parameter-set.rds")

#pyr_param_df_no_int <- pyr_param_df %>%
#  select(init_EIR, ivm_cov, Q0)
#
#saveRDS(pyr_param_df_no_int, "analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-no-int.rds")
##and save it to the DIDE drive.
#saveRDS(pyr_param_df_no_int, "W:/endectocides-cluster/data/scenario-parameter-set-no-int.rds")


#same but with different endec_mu and wane
#decide on grid space to explore

test_f <- function(endec,wane){
  t <- 1:365
  dd <- ifelse(t <= 23, round(0.132 + (endec * exp(-wane*t)),4), 0.132)
  test <- which(dd == 0.1320)[1] #the first day when value of dd comes back to 0.132
  return(test) #returns the day when dd first equals 0.132 or NA if it never happens
}
ff <- seq(0,1,0.01)
check <- array(dim = c(length(ff),length(ff)))
for(i in 1:length(ff)){ # i iterates over endec
  for(j in 1:length(ff)){ #j iterates over wane
    check[i,j] <- test_f(endec = ff[i], wane = ff[j]) #test f is called for each combination and stored in check[i,j].
    #check stores the days on which the function test_f() determines that the value of dd first returns to 0.132 for each combo of endec and wane parameters
  }
}
check #inspect check. Which array positions give a return to 30.
ff[2];ff[19] #this is the first one.

#what is minimum day when mu would return back to 0.132
require(tidyverse)
haz_curve <- read.table("IVM_derivation/ivermectin_hazards.txt", header = TRUE)

haz_curve <- haz_curve %>%
  select(day, d300)

ggplot(haz_curve, aes(x = day, y = d300))+
  geom_point()+
  geom_hline(aes(yintercept = 2), col = "red", lty = "dashed")+
  geom_hline(aes(yintercept = 1))+
  geom_vline(aes(xintercept = 23))+
  ylab("Hazard ratio")

which(haz_curve$d300 <= 2.00000)[1] #find first hazard ratio approaching 1

#first from day 15

#which values of i and j will give a return to mu = 0.132 FROM day 15
#find indices where condition is satisfied
#results <- which(check >= 15 & check <= 23, arr.ind = TRUE)
results <- which(check == 23, arr.ind = TRUE)
#results <- which(check == 24, arr.ind = TRUE) #which array position of check is equal to 30
results
#convert back to parameter values
endec_values <- ff[results[,1]]
wane_values <- ff[results[,2]]

output <- data.frame(endec = endec_values, wane = wane_values)
nrow(output)

range(output$endec) #0.01, 1
range(output$wane) #0.24, 0.70

nrow(output)
#index for the returning to d30
ff[2];ff[19] #0.01 and 0.18 - these are a bit off from what we think they should be
endec_mu_vec <- seq(0.01, 1, 0.11)
wane_vec <- seq(0.00,0.1, length.out = 10)

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

x <- pyr_param_df2_no_int %>%
  select(init_EIR, Q0, endec_mu, wane) %>%
  distinct()

saveRDS(pyr_param_df2_no_int, "analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-exp-decay-no-int.rds")
saveRDS(pyr_param_df2_no_int, "W:/endectocides-cluster/data/scenario-parameter-set-exp-decay-no-int.rds")


#pyr_param_df2 <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-exp-decay.rds")
#pyr_param_df2_no_int <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-exp-decay-no-int.rds")

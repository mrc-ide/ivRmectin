require(tidyverse)

#working out the correct endec_mu and wane values for malariasim_exp_decay.R

#read in odin models
odin_mda_early <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_early.rds")
odin_mda_medium <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_medium.rds")
odin_mda_late <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_late.rds")
odin_baseline <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_baseline.rds")
odin_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_ITN.rds")

#read in malariasim exp decay models
msim_mda_early <- readRDS("W:/endectocides-cluster/data/df_var1_mda_early_b_malariasim_odin.rds")
msim_mda_medium <- readRDS("W:/endectocides-cluster/data/df_var1_mda_medium_b_malariasim_odin.rds")
msim_mda_late <- readRDS("W:/endectocides-cluster/data/df_var1_mda_late_b_malariasim_odin.rds")
msim_baseline <- readRDS("W:/endectocides-cluster/data/df_baselineb_malariasim_odin.rds")
msim_ITN <- readRDS("W:/endectocides-cluster/data/df_df_var1_ITNb_malariasim_odin.rds")

#use least-squares fitting to calculate wane and endec_mu between models

#first filter in period of ivermectin MDA
time_period <- 365*15
mda_int <- 30

#introduce nets 1y into simulation

itn_on <- 365*5

net_seq <- seq(itn_on, time_period, by = 3*365)


#IVM MDA starts 6m after last ITN campaign
IVM_begin1 <- net_seq[2]+180
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

#IVM MDA starts 1y after last ITN campaign
IVM_begin2 <- net_seq[2] + 365
IVM_start2 <- c(IVM_begin2, IVM_begin2+mda_int, IVM_begin2+mda_int+mda_int)

#IVM starts 2.5 after last ITN campaign
IVM_begin3 <- net_seq[2] + (2*365)
IVM_start3 <- c(IVM_begin3, IVM_begin3+mda_int, IVM_begin3+mda_int+mda_int)

#parameter set for endec_mu and wane
endec_mu_vec <- seq(0,1,0.1)
wane_vec <- seq(0,1,0.1)
combos <- expand.grid(endec_mu_vec, wane_vec)


#early####
time_odin_mda_early <- odin_mda_early %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, Ivtot)

time_msim_mda_early <- msim_mda_early %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, Ivtot)

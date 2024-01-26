require(tidyverse)

#looking at how the endpoints (prevalence, EIR and mv) change across different bionomics settings
#in 2yo nets

#ivm is distributed at t = 3100, 3130 and 3160
#ivm killing for 23d so
ivm_on <- c(3100, 3130+23, 3160+23)
#look at the end points from ivm_on[1] to ivm_on[3]

#EXPLORE IN ANTAG MODEL FIRST

#first, where there is 10% resistance, by net type
antag_pyr <- read.csv("analysis/biting_LLINs/biting_antag_out_pyr_df.csv")
antag_pyr_pbo <- read.csv("analysis/biting_LLINs/biting_antag_out_pyr_pbo_df.csv")
antag_IG2 <- read.csv("analysis/biting_LLINs/biting_antag_out_IG2_df.csv")


no_ivm_pyr <- read.csv("analysis/biting_LLINs/pyr_out_df_antag_LLIN.csv", header = TRUE)

names(no_ivm_pyr)

####### looking in period of distribution only ######

antag_mods <- list.files(path = "C:/Users/nc1115/Documents/github/ivRmectin/analysis/biting_LLINs",
                         pattern = "antag", full.names = TRUE)


lapply(antag_mods, function(x){
  antag_ivm_on <- read.csv(x) %>%
    filter(between(t, ivm_on[1], ivm_on[3]))
  return(antag_ivm_on)
}) -> list_antag

#name them  by checking antag_mods
antag_mods
names(list_antag) <- c("df_antag_IG2", "df_antag_pyr", "df_antag_pyr_pbo")
list2env(list_antag, .GlobalEnv)
##################################################
dynamics_pyr <- antag_pyr %>%
  #filter(d_ITN0 == max(d_ITN0)) %>% # 10% resistance
  ggplot(aes(x = t, y =
               , col = as.factor(species)))+
  geom_line()+
  facet_grid(vars(itn_cov), vars(d_ITN0))
  #facet_wrap(~ itn_cov + d_ITN0)


#look at the time of the ivm distribution
no_ivm_pyr_df <- no_ivm_pyr %>%
  filter(itn_cov != 0.0 & between(t, ivm_on[1], ivm_on[3])) %>%
  group_by(itn_cov, d_ITN0, species) %>%
  summarise(mean_prev_LLIN = mean(slide_prev0to5),
            mean_EIR_LLIN = mean(EIR_tot),
            mean_mv_LLIN = mean(mv),
            mean_avhc_LLIN = mean(avhc),
            mean_mu_LLIN = mean(mu))

summary_antag_pyr <- antag_pyr %>%
  filter(itn_cov != 0.0 & between(t, ivm_on[1], ivm_on[3])) %>%
  group_by(itn_cov, d_ITN0, species) %>%
  summarise(mean_prev = mean(slide_prev0to5),
            mean_EIR = mean(EIR_tot),
            mean_mv = mean(mv),
            mean_avhc = mean(avhc),
            mean_mu = mean(mu))

summary2_antag_pyr <- left_join(summary_antag_pyr, no_ivm_pyr_df,  by = c("itn_cov", "d_ITN0", "species")) %>%
  mutate(red_prev = (mean_prev_LLIN - mean_prev)/mean_prev_LLIN,
         red_EIR = (mean_EIR_LLIN - mean_EIR)/mean_EIR_LLIN,
         red_mv = (mean_mv_LLIN - mean_mv)/mean_mv_LLIN,
         red_avhc = (mean_avhc_LLIN - mean_avhc)/mean_avhc_LLIN,
         red_mu = (mean_mu_LLIN - mean_mu)/mean_mu_LLIN)

prev_eff_pyr <- ggplot(summary2_antag_pyr, aes(x = itn_cov, y = red_prev, fill = as.factor(species)))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~d_ITN0)

eff_EIR_pyr <- ggplot(summary2_antag_pyr, aes(x = itn_cov, y = red_EIR, fill = as.factor(species)))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~d_ITN0)

eff_mv_pyr <- ggplot(summary2_antag_pyr, aes(x = itn_cov, y = red_mv, fill = as.factor(species)))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~d_ITN0)


species_names <- c("funestus" = "phi-B = 0.94",
                   "gambiae" = "phi-B = 0.85",
                   "arabiensis" = "phi-B = 0.71",
                   "stephensi" = "phi-B = 0.21")


time_between_df <- no_ivm_pyr_df
#FROM HERE: need to run again with the species for different bites_bed parameters

time_bt_meals <- time_between_df %>%
  group_by(species, itn_cov, d_ITN0) %>%
  summarise(time_bt_meals = 1/mean_avhc)


avhc_plot <- summary_antag_pyr |>
  mutate(across(species, ~factor(., levels = c("funestus", "gambiae", "arabiensis", "stephensi")))) %>%
  ggplot()+
  geom_bar(aes(x = itn_cov, y = 1/mean_avhc, fill = as.factor(d_ITN0)),
           stat = "identity", position = "dodge")+
  #facet_wrap(~species, labeller = as_labeller(species_names))+
  facet_wrap(~species)+
  #theme(legend.position = c(0.3, 0.85))+
  ylim(0, 5)+
  ylab("Average time between blood meals")


summary_antag_pyr |>
  mutate(across(species, ~factor(., levels = c("funestus", "gambiae", "arabiensis", "stephensi")))) %>%
  ggplot()+
  geom_bar(aes(x = species, y = 1/mean_avhc, fill = as.factor(itn_cov)),
           stat = "identity", position = "dodge")+
  facet_wrap(~d_ITN0)

avhc_plot_point <- summary_antag_pyr |>
  mutate(across(species, ~factor(., levels = c("funestus", "gambiae", "arabiensis", "stephensi")))) %>%
  ggplot()+
  geom_point(aes(x = itn_cov, y = 1/mean_avhc, col = as.factor(species)))+
  facet_wrap(~d_ITN0)

mu_plot <- summary_antag_pyr |>
  mutate(across(species, ~factor(., levels = c("funestus", "gambiae", "arabiensis", "stephensi")))) %>%
  ggplot()+
  geom_bar(aes(x = itn_cov, y = 1/mean_mu, fill = as.factor(d_ITN0)),
           stat = "identity", position = "dodge")+
  #facet_wrap(~species, labeller = as_labeller(species_names))+
  facet_wrap(~species)+
  guides(fill = "none")+
  ylab("Average life expectancy")

mu_plot_point <- summary_antag_pyr |>
  mutate(across(species, ~factor(., levels = c("funestus", "gambiae", "arabiensis", "stephensi")))) %>%
  ggplot()+
  geom_point(aes(x = itn_cov, y = 1/mean_mu, col = as.factor(species)))+
  facet_wrap(~d_ITN0)

ahvc_mu_plot <- cowplot::plot_grid(avhc_plot, mu_plot)

avhc_mu_plot_point <- cowplot::plot_grid(avhc_plot_point, mu_plot_point)

#phi-B moderately affects ivermectin uptake but more the mosquito life expectancy
#evening biting mosquitoes are as likely
#avhc is robust to resistance but mu is affected by it --> will affect mu-h


#why does red_EIR go negative? MODEL IS SAYING EIR HIGHER WITH NETS AND IVERMECTIN! maybe because lot of resistance? See ESA figure with neg EIR
#low coverage and low resistance --> by the time the net has aged,
summary2_antag_pyr %>%
  filter(itn_cov == 0.2) %>%
  select(mean_EIR_LLIN, mean_EIR, species)




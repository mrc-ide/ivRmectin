#plotting the dynamics of avhc and mu over net coverages and resistance levels

require(tidyverse)
require(cowplot)

breaks <- seq(2250, 3650, by = 30*6)
time_period <- 365*10
IVM_begin <- (365*8)+180
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)

eff_len <- 23
ivm_on <- IVM_start[1] #3100
ivm_off <- IVM_start[3]+eff_len #3183

antag_mod_LLIN_check <- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE) %>%
  filter(t == 3100) %>%
  group_by(itn_cov, d_ITN0, r_ITN0) %>%
  #summarise(mean_avhc = mean(avhc),
   #         mean_mu = mean(mu)) %>%
  mutate(s_ITN0 = 1 - d_ITN0 - r_ITN0)

#relationship between net efficacy and avhc looks odd.

antag_mod_LLIN<- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  group_by(itn_cov, d_ITN0, r_ITN0) %>%
  summarise(mean_avhc = mean(avhc),
           mean_mu = mean(mu)) %>%
  mutate(s_ITN0 = 1 - d_ITN0 - r_ITN0)

res_vec <- unique(antag_mod_LLIN$d_ITN0)

antag_mod_LLIN <- antag_mod_LLIN %>%
  mutate(res = case_when(d_ITN0 == res_vec[5] ~ "No resistance",
                         d_ITN0 == res_vec[4] ~ "10%",
                         d_ITN0 == res_vec[3] ~ "50%",
                         d_ITN0 == res_vec[2] ~ "70%",
                         d_ITN0 == res_vec[1] ~ "90%",
                         TRUE ~ NA_character_),
         plot = factor(res, levels = c("No resistance", "10%", "50%", "70%", "90%"),
                       labels =  c("No resistance", "10%", "50%", "70%", "90%"))) %>%
  filter(res != "No resistance", itn_cov %in% c(0.4, 0.6, 0.8))

pals <- c("#DDCC77", "#117733", "#332288", "#AA4499",
          "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

avhc_plot <- ggplot(antag_mod_LLIN, aes(x = itn_cov, y = 1/mean_avhc, fill = as.factor(res)))+
  geom_bar(stat = "identity", position = position_dodge())+
  ylim(0, 8)+
  theme_minimal()+
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8), minor_breaks = c(0.4, 0.6, 0.8), name = "Coverage (%)",
                     labels = c("40", "60", "80"))+
  scale_fill_manual(values = pals[4:1], name = "Resistance level")+
  ylab("Average time (days) between mosquito bloodmeals")+
  theme(legend.position = c(0.8, 0.9), legend.direction = "horizontal")

mu_plot <- ggplot(antag_mod_LLIN, aes(x = itn_cov, y = 1/mean_mu, fill = as.factor(res)))+
  geom_bar(stat = "identity", position = position_dodge())+
  ylim(0, 8)+
  theme_minimal()+
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8), minor_breaks = c(0.4, 0.6, 0.8), name = "Coverage (%)",
                     labels = c("40", "60", "80"))+
  scale_fill_manual(values = pals[4:1], name = "Resistance level")+
  ylab("Average mosquito life expectancy (days)")+
  theme(legend.position = "none")

biting_death_dynamics <- plot_grid(avhc_plot, mu_plot, labels = c("A)", "B)"))

ggsave(biting_death_dynamics, file = "plots/esa-figures/fig_4_biting_death.svg")


avhc_res_plot <- antag_mod_LLIN %>%
  filter(itn_cov == 0.8) %>%
  ggplot()+
  aes(x = plot, y = mean_avhc)+
  geom_bar(stat = "identity", position = position_dodge(), fill = "white", colour = "black")+
  theme_minimal()+
  #scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), minor_breaks = c(0, 0.2, 0.4, 0.6, 0.8), name = "Coverage (%)")+
  #scale_fill_manual(values = pals[5:1], name = "Resistance level")+
  ylab("Mean mosquito biting rate on humans")+
  ylim(0, 0.3)+
  xlab("Resistance")

mu_res_plot <- antag_mod_LLIN %>%
  filter(itn_cov == 0.8) %>%
  ggplot()+
  aes(x = plot, y = mean_mu)+
  geom_bar(stat = "identity", position = position_dodge(), fill = "white", colour = "black")+
  theme_minimal()+
  #scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), minor_breaks = c(0, 0.2, 0.4, 0.6, 0.8), name = "Coverage (%)")+
  #scale_fill_manual(values = pals[5:1], name = "Resistance level")+
  ylab("Mean mosquito mortality rate")+
  xlab("Resistance")+
  ylim(0, 0.3)

res_biting_dynamics <- plot_grid(avhc_res_plot, mu_res_plot, ncol = 2,
                                 labels = c("A)", "B)"))
ggsave(res_biting_dynamics, file = "plots/esa-figures/fig_4_res_biting_dynamics.svg")

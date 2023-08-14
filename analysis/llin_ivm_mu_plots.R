require(tidyverse)

#plotting fit of NC model to HS model ####

mv_fit_data<- read.csv(file = "data/llin_ivm_muh/mv_fit_data.csv",
                      header = TRUE)

time_period <- 365*10
IVM_begin <- (365*8)+180
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)

eff_len <- 23
ivm_on <- IVM_start[1] #3100
ivm_off <- IVM_start[3]+eff_len #3183



breaks_plot <- seq(2920, 3285, by = 30)

mv_plot_esa <- ggplot(mv_fit_data, aes(x = t, y = mv, col = as.factor(model)))+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "palegreen", alpha =0.25)+
  geom_line()+
  #xlim(2920, 3650)+
  ylab("Mosquito Density")+
  theme_minimal()+
  ylim(0, 50)+
  scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")+
  scale_colour_manual(values = c("black", "blue"), labels = c("Slater et al (2020) model: \n ivermectin mortality rate modelled with daily hazard ratios", "Modified Slater et al (2020) model: \n fixed additional mortality rate during ivermectin distribution period"),
                      name = "Modelling method")+
  theme(legend.position = c(.25, .95))+
  theme(text = element_text(size = 14))+
  #geom_vline(xintercept = c(ivm_on, ivm_off), linetype = "dashed")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 45, yend = 42, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 45, yend = 42, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 45, yend = 42, colour = "red", arrow = arrow())

ggsave(mv_plot_esa, file = "plots/llin_ivm_muh/mv_plot_esa.svg")


#reduction in EIR with different modelling methods and interventions####

#HS, IVM only
df_0_epi <- read.csv("data/llin_ivm_muh/df_0_epi.csv", header = TRUE)

#HS, IVM + NETS
df_3_epi <- read.csv("data/llin_ivm_muh/df_3_epi.csv", header = TRUE)


#HS, NETS only
df_3a_epi <- read.csv("data/llin_ivm_muh/df_3a_epi.csv", header = TRUE)

#HS, UTN and ivm
df_5_epi <- read.csv("data/llin_ivm_muh/df_5_epi.csv", header = TRUE)


#NC models

#NC, IVM only
df_2_epi <- read.csv("data/llin_ivm_muh/df_2_epi.csv", header = TRUE)


#NC, IVM + NETS
df_4_epi <- read.csv("data/llin_ivm_muh/df_4_epi.csv", header = TRUE)


#NC, NETS only
df_2a_epi <- read.csv("data/llin_ivm_muh/df_2a_epi.csv", header = TRUE)

#NC, IVM + UTN
df_6_epi <- read.csv("data/llin_ivm_muh/df_6_epi.csv", header = TRUE)

#absolute difference in EIR between df_0 and df_3
df_0_mean_EIR = df_0_epi %>% #ivm only
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 15.85293


df_3_mean_EIR = df_3_epi %>% #nets + ivm
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot))  #  12.2642 nets + ivm

df_5_mean_EIR = df_5_epi %>% #UTN + ivm
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot))  #  26.32482


df_0_mean_EIR - df_3_mean_EIR # 3.588733

df_0_mean_EIR - df_5_mean_EIR # -10.47189

df_2_mean_EIR = df_2_epi %>% #ivm only
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 14.32174 ivm only

df_4_mean_EIR = df_4_epi %>% # nets + ivm
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 11.54587. nets + ivm

df_2_mean_EIR - df_4_mean_EIR # 2.775873

df_6_mean_EIR <- df_6_epi %>% #untreated nets and ivm
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot)) #  23.22876

df_2_mean_EIR - df_6_mean_EIR #  -8.907015


#rel diff HS
rel_EIR_HS <- ((df_0_mean_EIR-df_3_mean_EIR)/df_0_mean_EIR)*100 #22.63766

#rel diff HS UTN
rel_EIR_UTN_HS <- ((df_0_mean_EIR-df_5_mean_EIR)/df_2_mean_EIR)*100 #-73.11886

#rel diff NC
rel_EIR_NC <- ((df_2_mean_EIR-df_4_mean_EIR)/df_2_mean_EIR)*100 #  19.38223

#rel diff NC untreated nets
rel_EIR_UTN_NC <- ((df_2_mean_EIR - df_6_mean_EIR)/df_2_mean_EIR)*100 #-62.19226

HS_plots_EIR <- ggplot(df_0_epi, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_line(data = df_3_epi, aes(x = t, y = EIR_tot), col = "red")+
  geom_line(data = df_3a_epi, aes(x = t, y = EIR_tot), col = "blue")+
  geom_line(data = df_5_epi, aes(x = t, y = EIR_tot), col = "purple")+
  ggtitle("EIR w HS model, ivm only (black), LLIN + IVM (red) and nets only (blue),
          untreated nets (purple)")

NC_plots_EIR <- ggplot(df_2_epi, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_line(data = df_4_epi, aes(x = t, y = EIR_tot), col = "red")+
  geom_line(data = df_2a_epi, aes(x = t, y = EIR_tot), col = "blue")+
  geom_line(data = df_6_epi, aes(x = t, y = EIR_tot), col = "purple")+
  ggtitle("EIR w NC model, ivm only (black), LLIN + IVM (red) and nets only (blue),
          untreated nets (purple)")

require(cowplot)
plot_grid(HS_plots_EIR, NC_plots_EIR)

df_0_epi_bind <- df_0_epi %>%
  mutate(model = "HS",
         int = "IVM only")

df_2_epi_bind <- df_2_epi %>%
  mutate(model = "NC",
         int = "IVM only")

df_4_epi_bind <- df_4_epi %>%
  mutate(model = "NC",
         int = "IVM + LLIN")

df_3_epi_bind <- df_3_epi %>%
  mutate(model = "HS",
         int = "IVM + LLIN")


df_5_epi_bind <- df_5_epi %>%
  mutate(model = "HS",
         int = "IVM only")

df_6_epi_bind <- df_6_epi %>%
  mutate(model = "NC",
         int = "IVM + LLIN")

ivm_llin_dat <- do.call("rbind",list(df_0_epi_bind,df_2_epi_bind,df_4_epi_bind, df_3_epi_bind))

ivm_utn_dat <- do.call("rbind", list(df_0_epi_bind, df_2_epi_bind, df_5_epi_bind, df_6_epi_bind))

ivm_llin_full_plot <- ggplot(ivm_llin_dat, aes(x = t, y = EIR_tot, col = as.factor(model), linetype = as.factor(int)))+
  geom_line()+
  theme_minimal()+
 # scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")+
  ylab("EIR")+
  labs(linetype = "Intervention")+
  scale_colour_manual(values = c("black", "blue"), labels = c("Slater et al (2020) model: \n ivermectin mortality rate modelled with daily hazard ratios and intervention-dependent feeding rate", "Modified Slater et al (2020) model: \n fixed additional mortality rate during ivermectin distribution period and fixed feeding rate"),
                      name = "Modelling method")+
  theme(legend.position = c(0.5, 0.8), text = element_text(size = 14))+
  xlab("Time (days)")

ivm_llin_red_plot <- ggplot(ivm_llin_dat, aes(x = t, y = EIR_tot, col = as.factor(model), linetype = as.factor(int)))+
  geom_line()+
  theme_minimal()+
  scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")+
  ylab("EIR")+
  #labs(linetype = "Intervention")+
  scale_colour_manual(values = c("black", "blue"), labels = c("Slater et al (2020) model: \n ivermectin mortality rate modelled with daily hazard ratios and intervention-dependent feeding rate", "Modified Slater et al (2020) model: \n fixed additional mortality rate during ivermectin distribution period and fixed feeding rate"),
                      name = "Modelling method")+
  theme(legend.position = "none", text = element_text(size = 14))

ivm_llin_eir_plot <- plot_grid(ivm_llin_full_plot, ivm_llin_red_plot, nrow = 2)

ggsave(ivm_llin_eir_plot,
       file = "plots/llin_ivm_muh/ivm_llin_eir_plot.svg")

#untreated nets
ivm_utn_full_plot <-ggplot(ivm_llin_dat, aes(x = t, y = EIR_tot, col = as.factor(model), linetype = as.factor(int)))+
  geom_line()+
  theme_minimal()+
  # scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")+
  ylab("EIR")+
  #labs(linetype = "Intervention")+
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("IVM only", "IVM + UTN"), name = "Intervention")+
  scale_colour_manual(values = c("black", "blue"), labels = c("Slater et al (2020) model: \n ivermectin mortality rate modelled with daily hazard ratios and intervention-dependent feeding rate", "Modified Slater et al (2020) model: \n fixed additional mortality rate during ivermectin distribution period and fixed feeding rate"),
                      name = "Modelling method")+
  theme(legend.position = c(0.5, 0.8), text = element_text(size = 14))+
  xlab("Time (days)")

ivm_utn_red_plot <- ggplot(ivm_llin_dat, aes(x = t, y = EIR_tot, col = as.factor(model), linetype = as.factor(int)))+
  geom_line()+
  theme_minimal()+
  ylab("EIR")+
  #labs(linetype = "Intervention")+
  #scale_linetype_manual(values = c("solid", "dashed"), labels = c("IVM only", "IVM + UTN"), name = "Intervention")+
  #scale_colour_manual(values = c("black", "blue"), labels = c("Slater et al (2020) model: \n ivermectin mortality rate modelled with daily hazard ratios and intervention-dependent feeding rate", "Modified Slater et al (2020) model: \n fixed additional mortality rate during ivermectin distribution period and fixed feeding rate"),
  #                    name = "Modelling method")+
  theme(legend.position = c(0.5, 0.8), text = element_text(size = 14))+
  #xlab("Time (days)")
  scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")

ivm_utn_eir_plot <- plot_grid(ivm_utn_full_plot, ivm_utn_red_plot, nrow = 2)

ivm_eir_plots <- plot_grid(ivm_llin_full_plot, ivm_utn_full_plot, ivm_llin_red_plot,  ivm_utn_red_plot)
ggsave(ivm_eir_plots, file = "plots/llin_ivm_muh/ivm_eir_plots.svg")

#proportion of mosquitoes killed by ivermectin under different LLIN coverages and endophagy levels####

HS_out_df <- read.csv(file = "data/llin_ivm_muh/HS_itn_cov_loop.csv", header = TRUE)
NC_out_df <- read.csv(file = "data/llin_ivm_muh/NC_itn_cov_loop.csv", header = TRUE)

HS_out_df_utn <- read.csv(file = "data/llin_ivm_muh/HS_utn_cov_loop.csv", header = TRUE)
NC_out_df_utn <- read.csv(file = "data/llin_ivm_muh/NC_utn_cov_loop.csv", header = TRUE)

#across itn coverages in the HS model, how does the proportion of ivm-killed mosquitoes change?

HS_prop_killed_ivm_summary <- HS_out_df %>%
  group_by(itn_cov) %>%
  summarise(mv_dead_tot = sum(mv_dead),
            mvx_dead_tot = sum(mvx_dead)) %>%
  mutate(prop_killed_ivm_tot = (mvx_dead_tot/mv_dead_tot)*100) %>%
  mutate(rel_diff_prop_killed_ivm = (prop_killed_ivm_tot[1] - prop_killed_ivm_tot)/prop_killed_ivm_tot[1],
         model = "HS",
         int = "LLIN")


NC_prop_killed_ivm_summary <- NC_out_df %>%
  group_by(itn_cov) %>%
  summarise(mv_dead_tot = sum(mv_dead),
            mvx_dead_tot = sum(mvx_dead)) %>%
  mutate(prop_killed_ivm_tot = (mvx_dead_tot/mv_dead_tot)*100) %>%
  mutate(rel_diff_prop_killed_ivm = (prop_killed_ivm_tot[1] - prop_killed_ivm_tot)/prop_killed_ivm_tot[1],
         model ="NC",
         int = "LLIN")


prop_killed_ivm_summary_llin <- rbind(HS_prop_killed_ivm_summary, NC_prop_killed_ivm_summary)

prop_killed_ivm_plot_llin <- ggplot(prop_killed_ivm_summary_llin, aes(x = itn_cov, y = prop_killed_ivm_tot, fill = as.factor(model)))+
  geom_bar(stat = "identity", position = position_dodge())+
  ylim(0, 1)+
  theme_minimal()+
  labs(x = "LLIN coverage (%)", y = "Proportion of mosquitoes \n killed by ivermectin",
       fill = "Model type")

#repeat for utn

HS_prop_killed_ivm_summary_utn <- HS_out_df_utn %>%
  group_by(itn_cov) %>%
  summarise(mv_dead_tot = sum(mv_dead),
            mvx_dead_tot = sum(mvx_dead)) %>%
  mutate(prop_killed_ivm_tot = (mvx_dead_tot/mv_dead_tot)*100) %>%
  mutate(rel_diff_prop_killed_ivm = (prop_killed_ivm_tot[1] - prop_killed_ivm_tot)/prop_killed_ivm_tot[1],
         model = "HS",
         int = "UTN")


NC_prop_killed_ivm_summary_utn <- NC_out_df_utn %>%
  group_by(itn_cov) %>%
  summarise(mv_dead_tot = sum(mv_dead),
            mvx_dead_tot = sum(mvx_dead)) %>%
  mutate(prop_killed_ivm_tot = (mvx_dead_tot/mv_dead_tot)*100) %>%
  mutate(rel_diff_prop_killed_ivm = (prop_killed_ivm_tot[1] - prop_killed_ivm_tot)/prop_killed_ivm_tot[1],
         model ="NC",
         int = "UTN")


prop_killed_ivm_summary_utn <- rbind(HS_prop_killed_ivm_summary_utn, NC_prop_killed_ivm_summary_utn)

prop_killed_ivm_plot_utn <- ggplot(prop_killed_ivm_summary_utn, aes(x = itn_cov, y = prop_killed_ivm_tot, fill = as.factor(model)))+
  geom_bar(stat = "identity", position = position_dodge())+
  ylim(0, 1)+
  theme_minimal()+
  labs(x = "UTN coverage (%)", y = "Proportion of mosquitoes \n killed by ivermectin",
       fill = "Model type")

prop_killed_net_plots <- plot_grid(prop_killed_ivm_plot_llin, prop_killed_ivm_plot_utn, ncol = 2)

prop_killed_all <- do.call("rbind", list(HS_prop_killed_ivm_summary, HS_prop_killed_ivm_summary_utn,
                                         NC_prop_killed_ivm_summary, NC_prop_killed_ivm_summary_utn))

prop_killed_all <- prop_killed_all %>%
  mutate(model_int = paste(model, int))

prop_killed_ivm_plot <- ggplot(prop_killed_all, aes(x = itn_cov, y = prop_killed_ivm_tot , fill = as.factor(model_int)))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(name = "Model and intervention type", values = c("black", "grey",
                                                                     "blue", "turquoise"),
                    labels = c("HS model, LLIN", "HS model, UTN", "NC moodel, LLIN", "NC model, UTN"))+
  theme_minimal()+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), minor_breaks = c(0, 0.2, 0.4, 0.6, 0.8), name = "Coverage (%)")+
  ylab("Proportion of mosquites killed by ivermectin")+
  ylim(0, 0.5)+
  theme(legend.position = c(0.85, 0.8))
ggsave(prop_killed_ivm_plot, file = "plots/llin_ivm_muh/prop_killed_ivm_plot.svg")

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

HS_ivm_only <- read.csv("data/llin_ivm_muh/df_0_epi.csv")
NC_ivm_only <- read.csv("data/llin_ivm_muh/df_2_epi.csv")


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


#plot showing EIR and reduction with IVM only vs IVM + LLIN
pyr_HS_df <- read.csv("data/llin_ivm_muh/pyr_out_df_HS.csv", header = TRUE) %>%
  filter(itn_cov == 0.8)


breaks_plot <- seq(2920, 3285, by = 30)
res_vector <- c(0, 0.1, 0.5, 0.7, 0.9)
cbPalette <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73", "#999999", "#CC79A7")
cbpPalette2 <- c("black", cbPalette[1:5])

HS_ivm_LLIN_dynamics <- ggplot(HS_ivm_only, aes(x = t, y = EIR_tot))+
  geom_line(aes(colour = "black"))+
  geom_line(data = pyr_HS_df, aes(x = t, y = EIR_tot, colour = as.factor(d_ITN0)))+
  ylim(-1, 100)+
  scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = -1, yend = 4, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = -1, yend = 4, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = -1, yend = 4, colour = "red", arrow = arrow())+
  theme_minimal()+
  ylab("EIR (HS model)")+
  scale_colour_manual(values = c(cbpPalette2), labels = c("IVM only", "IVM + LLIN (0% resistance)",
                                                          "IVM + LLIN (10% resistance)",
                                                          "IVM + LLIN (50% resistance)",
                                                          "IVM + LLIN (70% resistance)",
                                                          "IVM + LLIN (90% resistance)"),
                      name = "Intervention")+
  theme(legend.position = "none")

#ggsave(HS_ivm_LLIN_dynamics, file = "plots/llin_ivm_muh/HS_ivm_llin_dynamics_plot.svg")

pyr_NC_df <- read.csv("data/llin_ivm_muh/pyr_out_df_NC.csv", header = TRUE) %>%
  filter(itn_cov == 0.8)

NC_ivm_LLIN_dynamics <- ggplot(NC_ivm_only, aes(x = t, y = EIR_tot))+
  geom_line(aes(colour = "black"))+
  geom_line(data = pyr_NC_df, aes(x = t, y = EIR_tot, colour = as.factor(d_ITN0)))+
  ylim(-1, 100)+
  scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = -1, yend = 4, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = -1, yend = 4, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = -1, yend = 4, colour = "red", arrow = arrow())+
  theme_minimal()+
  ylab("EIR (NC model)")+
  scale_colour_manual(values = c(cbpPalette2), labels = c("IVM only", "IVM + LLIN (0% resistance)",
                                                          "IVM + LLIN (10% resistance)",
                                                          "IVM + LLIN (50% resistance)",
                                                          "IVM + LLIN (70% resistance)",
                                                          "IVM + LLIN (90% resistance)"),
                      name = "Intervention")+
  theme(legend.position = c(0.7, 0.9))
require(cowplot)
ivm_llin_dynamics <- plot_grid(HS_ivm_LLIN_dynamics, NC_ivm_LLIN_dynamics, labels = c("A", "B"))
ggsave(ivm_llin_dynamics, file = "plots/llin_ivm_muh/ivm_llin_dynamics.svg")


#proportion of dead mosquitoes across different net coverages
#compare HS and NC model


prop_dead_HS <- read.csv("data/llin_ivm_muh/pyr_out_df_HS.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off), d_ITN0 > 0.34) %>% #no resistance
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead),
            model = "HS_pyr")

prop_dead_NC <- read.csv("data/llin_ivm_muh/pyr_out_df_NC.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off), d_ITN0 > 0.34) %>% #no resistance
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead),
            model = "NC_pyr")


#illustrate proportion of mosquitoes dying over time with HS model
prop_killed_ivm_hs <- ggplot(prop_dead_HS, aes(x = itn_cov, y = prop_dead_ivm))+
  geom_bar(stat = "identity", position = position_dodge(), fill = "black")+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), minor_breaks = c(0, 0.2, 0.4, 0.6, 0.8), name = "LLIN Coverage (%)")+
  ylab("Proportion of mosquites killed by ivermectin")+
  theme_minimal()+
  ylim(0, 0.01)
ggsave(prop_killed_ivm_hs, file = "plots/llin_ivm_muh/prop_killed_ivm_hs.svg")

utn_df_HS <- read.csv("data/llin_ivm_muh/utn_out_df_HS.csv") %>%
  group_by(itn_cov) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "HS_utn")


utn_df_NC <- read.csv("data/llin_ivm_muh/utn_out_df_NC.csv") %>%
  group_by(itn_cov) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "NC_utn")

utn_df <- rbind(utn_df_HS, utn_df_NC)


prop_killed_all <- rbind(prop_dead_HS, prop_dead_NC, utn_df)

prop_killed_int <- ggplot(prop_killed_all, aes(x = itn_cov, y = prop_dead_ivm , fill = as.factor(model)))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(name = "Model and intervention type", values = c("black", "grey",
                                                                     "blue", "turquoise"))+
  theme_minimal()+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), minor_breaks = c(0, 0.2, 0.4, 0.6, 0.8), name = "Coverage (%)")+
  ylab("Proportion of mosquites killed by ivermectin")+
  ylim(0, 0.01)+
  theme(legend.position = c(0.65, 0.8))

ggsave(prop_killed_int, file = "plots/llin_ivm_muh/prop_killed_int.svg")

ggplot(utn_df, aes(x = itn_cov, y = prop_dead_ivm, fill = as.factor(model)))+
  geom_bar(stat = "identity", position = position_dodge())+
  ylim(0, 0.01)

#unique(pyr_df_HS$ref)
#unique(pyr_df_HS$d_ITN0)
pyr_df_HS <- read.csv("data/llin_ivm_muh/pyr_out_df_HS.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "HS pyr")

pyr_df_NC <- read.csv("data/llin_ivm_muh/pyr_out_df_NC.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "NC pyr")

pyr_df <- rbind(pyr_df_HS, pyr_df_NC)

summary_ivm_HS <- HS_ivm_only %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead))

summary_ivm_NC <- NC_ivm_only %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead))


pyr_df_HS_plot <- ggplot(pyr_df_HS, aes(x = itn_cov, y = mean_EIR, fill = as.factor(d_ITN0)))+
  geom_bar(stat = "identity", position = position_dodge())

pyr_df_NC_plot <- ggplot(pyr_df_NC, aes(x = itn_cov, y = mean_EIR, fill = as.factor(d_ITN0)))+
  geom_bar(stat = "identity", position = position_dodge())

pyr_summary <- pyr_df %>%
  filter(d_ITN0 > 0.34) %>%
  mutate(base_ivm = case_when(grepl("NC", model) ~ summary_ivm_NC$mean_EIR,
                              grepl("HS", model, ignore.case = TRUE) ~ summary_ivm_HS$mean_EIR),
         abs_diff_ivm = mean_EIR - base_ivm,
         rel_diff_ivm = (mean_EIR - base_ivm)/mean_EIR)

utn_summary <- as.data.frame(utn_df) %>%
  mutate(d_ITN0 = rep(0.059, times = 10),
         base_ivm = case_when(grepl("NC", model) ~ summary_ivm_NC$mean_EIR,
                              grepl("HS", model, ignore.case = TRUE) ~ summary_ivm_HS$mean_EIR),
         abs_diff_ivm = mean_EIR - base_ivm,
         rel_diff_ivm = (mean_EIR - base_ivm)/mean_EIR,
         res = "utn")

utn_summary <- utn_summary[, c("itn_cov", "d_ITN0", "mean_EIR", "prop_dead_ivm", "model", "base_ivm",
                               "abs_diff_ivm", "rel_diff_ivm", "res")]

res_val <- rep(c(0, 0.1, 0.5, 0.7, 0.9), times = 10)

pyr_summary_res <- pyr_df %>%
  #filter(d_ITN0 > 0.34) %>%
  mutate(base_ivm = case_when(grepl("NC", model) ~ summary_ivm_NC$mean_EIR,
                              grepl("HS", model, ignore.case = TRUE) ~ summary_ivm_HS$mean_EIR),
         abs_diff_ivm = mean_EIR - base_ivm,
         rel_diff_ivm = (mean_EIR - base_ivm)/mean_EIR) %>%
  ungroup () %>%
  mutate(res = res_val)

IG2_df_HS <- read.csv("data/llin_ivm_muh/IG2_out_df_HS.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "HS IG2")

IG2_df_NC <- read.csv("data/llin_ivm_muh/IG2_out_df_NC.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "NC IG2")

IG2_df <- rbind(IG2_df_HS, IG2_df_NC)

IG2_summary_res <- IG2_df %>%
  mutate(base_ivm = case_when(grepl("NC", model) ~ summary_ivm_NC$mean_EIR,
                              grepl("HS", model, ignore.case = TRUE) ~ summary_ivm_HS$mean_EIR),
         abs_diff_ivm = mean_EIR - base_ivm,
         rel_diff_ivm = (mean_EIR - base_ivm)/mean_EIR) %>%
  ungroup () %>%
  mutate(res = res_val)

pbo_df_HS <- read.csv("data/llin_ivm_muh/pyr_pbo_out_df_HS.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "HS PBO")

pbo_df_NC <- read.csv("data/llin_ivm_muh/pyr_pbo_out_df_NC.csv") %>%
  group_by(itn_cov, d_ITN0) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            prop_dead_ivm = sum(mvx_dead)/sum(mv_dead)) %>%
  mutate(model = "NC PBO")

pbo_df <- rbind(pbo_df_HS, pbo_df_NC)


pbo_summary_res <- pbo_df %>%
  mutate(base_ivm = case_when(grepl("NC", model) ~ summary_ivm_NC$mean_EIR,
                              grepl("HS", model, ignore.case = TRUE) ~ summary_ivm_HS$mean_EIR),
         abs_diff_ivm = mean_EIR - base_ivm,
         rel_diff_ivm = (mean_EIR - base_ivm)/mean_EIR) %>%
  ungroup () %>%
  mutate(res = res_val)

net_summary <- do.call("rbind", list(pyr_summary_res, utn_summary,  pbo_summary_res, IG2_summary_res))

res_vector <- c(0, 0.1, 0.5, 0.7, 0.9)
#cbPalette <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
bioassay_mort <- 1-res_vector

net_summary_facet <- net_summary %>%
  separate(model, c("model", "net"))

names <- list(
  "HS" = "HS model",
  "NC" = "NC model",
  "IG2" = "IG2 net",
  "PBO" = "PBO net",
  "pyr" = "Pyrethroid-only net",
  "utn" = "UTN"
)

net_labeller <- function(variable, value){
  return(names[value])
}

#across coverages and resistance profiles
red_eir_nets <- ggplot(net_summary_facet, aes(x = itn_cov, y = rel_diff_ivm*100, fill = res))+
  geom_bar(stat = "identity", position = position_dodge())+
  #facet_grid(~fct_relevel(model,"HS pyr", "NC pyr", "HS PBO", "NC PBO", "HS IG2", "NC IG2", "HS UTN", "NC UTN"), cols = 2, rows = 3)+
 # xlab("Net coverage (%)")+
  facet_grid(cols = vars(model), rows = vars(net),
             labeller = net_labeller)+
  #facet_wrap(~fct_relevel())
  ylab("Relative difference in EIR (%)")+
  scale_fill_manual(name = "Resistance", values = cbPalette, labels = c("0%", "10%", "50%", "70%", "90%", "UTN"))+
  #scale_fill_viridis_b(name = "Resistance", labels = c("0", "10%", "50%", "70%", "90%", "UTN") )+
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), minor_breaks = c(0, 0.2, 0.4, 0.6, 0.8), name = "Coverage (%)")+
  theme_minimal()+
  theme(panel.spacing = unit(2, "lines"))+
  geom_hline(yintercept = 0, linetype = "solid", colour = "black")+
  theme(legend.position = c(0.7, 0.1), legend.direction = "horizontal")
ggsave(red_eir_nets, file = "plots/llin_ivm_muh/red_eir_nets.svg")

#avhc and mu plots

pyr_HS <- read.csv("data/llin_ivm_muh/pyr_out_df_HS.csv")  %>%
  mutate(model = "HS")
pyr_NC <- read.csv("data/llin_ivm_muh/pyr_out_df_NC.csv") %>%
  mutate(model = "NC")

res_pyr <- unique(pyr_HS$d_ITN0)

pyr_HS <- pyr_HS %>%
  mutate(res = case_when(d_ITN0 == res_pyr[1] ~ 0,
                         d_ITN0 == res_pyr[2] ~ 0.1,
                         d_ITN0 == res_pyr[3] ~ 0.5,
                         d_ITN0 == res_pyr[5] ~ 0.7,
                         d_ITN0 == res_pyr[7] ~ 0.9,
                         TRUE ~ NA_real_))

pyr_NC <- pyr_NC %>%
  mutate(res = case_when(d_ITN0 == res_pyr[1] ~ 0,
                         d_ITN0 == res_pyr[2] ~ 0.1,
                         d_ITN0 == res_pyr[3] ~ 0.5,
                         d_ITN0 == res_pyr[5] ~ 0.7,
                         d_ITN0 == res_pyr[7] ~ 0.9,
                         TRUE ~ NA_real_))

utn_HS <- read.csv("data/llin_ivm_muh/utn_out_df_HS.csv") %>%
  mutate(model = "HS")
utn_NC <- read.csv("data/llin_ivm_muh/utn_out_df_NC.csv")%>%
  mutate(model = "NC")


utn_HS <- utn_HS %>%
  mutate(res = 99)

utn_NC <- utn_NC %>%
  mutate(res = 99)


pyr_pbo_HS <- read.csv("data/llin_ivm_muh/pyr_pbo_out_df_HS.csv") %>%
  mutate(model = "HS")
pyr_pbo_NC <- read.csv("data/llin_ivm_muh/pyr_pbo_out_df_NC.csv")%>%
  mutate(model = "NC")

res_pbo <- unique(pyr_pbo_HS$d_ITN0)

pyr_pbo_HS <- pyr_pbo_HS %>%
  mutate(res = case_when(d_ITN0 == res_pbo[1] ~ 0,
                         d_ITN0 == res_pbo[2] ~ 0.1,
                         d_ITN0 == res_pbo[3] ~ 0.5,
                         d_ITN0 == res_pbo[5] ~ 0.7,
                         d_ITN0 == res_pbo[7] ~ 0.9,
                         TRUE ~ NA_real_))
pyr_pbo_NC <- pyr_pbo_NC %>%
  mutate(res = case_when(d_ITN0 == res_pbo[1] ~ 0,
                         d_ITN0 == res_pbo[2] ~ 0.1,
                         d_ITN0 == res_pbo[3] ~ 0.5,
                         d_ITN0 == res_pbo[5] ~ 0.7,
                         d_ITN0 == res_pbo[7] ~ 0.9,
                         TRUE ~ NA_real_))

IG2_HS <- read.csv("data/llin_ivm_muh/IG2_out_df_HS.csv") %>%
  mutate(model = "HS")
IG2_NC <- read.csv("data/llin_ivm_muh/IG2_out_df_NC.csv")%>%
  mutate(model = "NC")

res_ig2 <- unique(IG2_HS$d_ITN0)

IG2_HS <- IG2_HS %>%
  mutate(res = case_when(d_ITN0 == res_ig2[1] ~ 0,
                         d_ITN0 == res_ig2[2] ~ 0.1,
                         d_ITN0 == res_ig2[3] ~ 0.5,
                         d_ITN0 == res_ig2[5] ~ 0.7,
                         d_ITN0 == res_ig2[7] ~ 0.9,
                         TRUE ~ NA_real_))

IG2_NC <- IG2_NC %>%
  mutate(res = case_when(d_ITN0 == res_ig2[1] ~ 0,
                         d_ITN0 == res_ig2[2] ~ 0.1,
                         d_ITN0 == res_ig2[3] ~ 0.5,
                         d_ITN0 == res_ig2[5] ~ 0.7,
                         d_ITN0 == res_ig2[7] ~ 0.9,
                         TRUE ~ NA_real_))


nets <- do.call("rbind", list(pyr_HS, pyr_NC, utn_HS, utn_NC,
                              pyr_pbo_HS, pyr_pbo_NC, IG2_HS, IG2_NC))

net_biting <- nets %>%
  group_by(model, net_type, itn_cov, res) %>%
  summarise(mean_avhc = mean(avhc),
            mean_mu = mean(mu),
            res = res) %>%
  ggplot()+
  aes(x = itn_cov, y = mean_avhc, colour = as.factor(res)) +
  geom_point()+
  geom_point(aes(x = itn_cov, y = mean_mu, colour = as.factor(res)), shape = 4)+
  facet_grid(cols = vars(model), rows = vars(net_type))+
  scale_fill_manual(name = "Resistance", values = cbPalette, labels = c("0%", "10%", "50%", "70%", "90%", "UTN"))


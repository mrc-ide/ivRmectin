require(tidyverse)
require(cowplot)

#remove these times because there are EIR spikes here
rm <- c(2290, 2291, 2292, 3385, 3386, 3387)

antag_mod_LLIN<- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov) %>%
  mutate(model = "antagonistic", int = "LLIN") %>%
  filter(!t %in% rm)

add_mod_LLIN<- read.csv("analysis/esa-analysis/output/pyr_out_add.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov) %>%
  mutate(model = "additive", int = "LLIN") %>%
  filter(!t %in% rm)

antag_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_antag_IVM.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov) %>%
  mutate(model = "antagonistic", int = "LLIN & IVM") %>%
  filter(!t %in% rm)

add_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_add_IVM.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov) %>%
  mutate(model = "additive", int = "LLIN & IVM") %>%
  filter(!t %in% rm)


pals <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
          "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#breaks_plot <- seq(2250, 3650, by = 30)
breaks <- seq(2250, 3650, by = 30*6)
time_period <- 365*10
IVM_begin <- (365*8)+180
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)

eff_len <- 23
ivm_on <- IVM_start[1] #3100
ivm_off <- IVM_start[3]+eff_len #3183

net_seq <- seq(100, 3650, by = 3*365)

death_params <- unique(antag_mod_LLIN$d_ITN0)

antag_mod <- rbind(antag_mod_LLIN, antag_mod_LLIN_IVM) %>%
  mutate(res = case_when(d_ITN0 == death_params[1] ~ "No resistance",
                         d_ITN0 == death_params[2] ~ "10% resistance",
                         d_ITN0 == death_params[3] ~ "50% resistance",
                         d_ITN0 == death_params[4] ~ "70% resistance",
                         d_ITN0 == death_params[5] ~ "90% resistance",
                         TRUE ~ NA_character_)) %>%
  filter(res != "No resistance")



add_mod <- rbind(add_mod_LLIN, add_mod_LLIN_IVM) %>%
  mutate(res = case_when(d_ITN0 == death_params[1] ~ "No resistance",
                         d_ITN0 == death_params[2] ~ "10% resistance",
                         d_ITN0 == death_params[3] ~ "50% resistance",
                         d_ITN0 == death_params[4] ~ "70% resistance",
                         d_ITN0 == death_params[5] ~ "90% resistance",
                         TRUE ~ NA_character_)) %>%
  filter(res != "No resistance")

pals2 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
          "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#filter to the no resistance plot
mv_80_cov_no_res <- antag_mod %>%
  filter(itn_cov == 0.8 & res == "10% resistance") %>%
  ggplot()+
  aes(x = t, y = mv, colour = as.factor(int))+
  geom_line()+
  scale_colour_manual(values = c(pals2[2], pals2[1]), name = "Intervention")+
  theme_minimal()+
  ylab("Mosquito Density \n (relative to human population size)")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 45, yend = 38, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 45, yend = 38, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 45, yend = 38, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 45, yend = 38, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 45, yend = 38, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  theme(legend.position = c(0.8, 0.3))


mv_llin_antag <- antag_mod %>%
  filter(itn_cov == 0.8 & between(t, ivm_on, ivm_off)) %>%
  group_by(res, int) %>%
  summarise(mean_mv = mean(mv))

mv_llin_antag_wide <- mv_llin_antag %>%
  spread(int, mean_mv)

mv_llin_antag_summary <- mv_llin_antag_wide %>%
  mutate(abs_diff = LLIN - `LLIN & IVM`,
         rel_diff = ((LLIN - `LLIN & IVM`)/LLIN)*100)

res_order <-  c("10% resistance","50% resistance","70% resistance","90% resistance")

mv_llin_add <- add_mod %>%
  filter(itn_cov == 0.8 & between(t, ivm_on, ivm_off)) %>%
  group_by(res, int) %>%
  summarise(mean_mv = mean(mv))

mv_llin_add_wide <- mv_llin_add %>%
  spread(int, mean_mv)

mv_llin_add_summary <- mv_llin_add_wide %>%
  mutate(abs_diff = LLIN - `LLIN & IVM`,
         rel_diff = ((LLIN - `LLIN & IVM`)/LLIN)*100)

#plot the relative differences
rel_diff_antag <- ggplot(mv_llin_antag_summary, aes(x = res, y = rel_diff))+
  geom_bar(stat = "identity", position = position_dodge(), colour = "black",fill = "white")+
  ylim(0, 30)+
  theme_minimal()+
  labs(x = "Resistance", y = "Relative difference \n in mosquito density (%) \n due to ivermectin")+
  scale_x_discrete(limits = res_order)+
  ggtitle("Antagonistic model")

rel_diff_add <- ggplot(mv_llin_add_summary, aes(x = res, y = rel_diff))+
  geom_bar(stat = "identity", position = position_dodge(), colour = "black",fill = "white")+
  ylim(0, 30)+
  theme_minimal()+
  labs(x = "Resistance", y = "Relative difference \n in mosquito density (%) \n due to ivermectin")+
  scale_x_discrete(limits = res_order)+
  ggtitle("Additive model")

rel_diff_mv <- plot_grid(rel_diff_antag, rel_diff_add, labels = c("B)", "C)"), label_x = 0.05)

density_plots <- plot_grid(mv_80_cov_no_res, rel_diff_mv, labels = c("A)", ""), nrow =2,
                           label_x = 0.05)
ggsave(density_plots, file = "plots/esa-figures/fig_2_density_plots.svg")

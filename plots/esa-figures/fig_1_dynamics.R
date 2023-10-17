#dynamics plot

#plot prevalence and EIR for both models

require(tidyverse)
require(cowplot)

df <- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE)
res <- unique(df$d_ITN0)



spikes <-  read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE)

rm1 <- spikes %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2]) %>%
  filter(between(t, 1190, 1199)) #remove t = 1195, t = 1196, t = 1197

spikes %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2]) %>%
  filter(between(t, 3375, 3395)) #remove t = 2290, t = 2291, t = 2292

rm <- c(2290, 2291, 2292, 3385, 3386, 3387)

#read in the other files
antag_mod_LLIN <- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE) %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov, EIR_tot, slide_prev0to5) %>%
  mutate(model = "antagonistic", int = "LLIN", id = row_number())

add_mod_LLIN<- read.csv("analysis/esa-analysis/output/pyr_out_add.csv", header = TRUE) %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov,EIR_tot, slide_prev0to5) %>%
  mutate(model = "additive", int = "LLIN", id = row_number())

antag_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_antag_IVM.csv", header = TRUE) %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov,EIR_tot, slide_prev0to5) %>%
  mutate(model = "antagonistic", int = "LLIN & IVM", id = row_number())

add_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_add_IVM.csv", header = TRUE) %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov,EIR_tot, slide_prev0to5) %>%
  mutate(model = "additive", int = "LLIN & IVM", id = row_number())

add_mod <- rbind(add_mod_LLIN, add_mod_LLIN_IVM)
antag_mod <- rbind(antag_mod_LLIN, antag_mod_LLIN_IVM)

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

add_EIR <- ggplot(add_mod, aes(x = t, y = EIR_tot, colour = as.factor(int)))+
  geom_line()+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  scale_colour_manual(values = c(pals[2], pals[1]), name = "Intervention") +
  theme(legend.position = "none", text = element_text(size = 12))+
  ylab("EIR")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 50, yend = 45, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 50, yend = 45, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 50, yend = 45, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 75, yend = 65, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 75, yend = 65, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 80)+
  ggtitle("Additive model")


antag_EIR <- ggplot(antag_mod, aes(x = t, y = EIR_tot, colour = as.factor(int)))+
  geom_line()+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  scale_colour_manual(values = c(pals[2], pals[1]), name = "Intervention") +
  theme(legend.position  = "none", text = element_text(size = 12))+
  ylab("EIR")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 50, yend = 45, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 50, yend = 45, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 50, yend = 45, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 75, yend = 65, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 75, yend = 65, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 80)+
  ggtitle("Antagonistic model")



EIR_plot <- plot_grid(antag_EIR, add_EIR, ncol = 2,
                      labels = c("A)", "B)"))




#ggsave(EIR_plot, file = "plots/esa-figures/fig_1_dynamics.svg")


#prevalence
add_prev <- ggplot(add_mod, aes(x = t, y = slide_prev0to5*100, colour = as.factor(int)))+
  geom_line()+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  scale_colour_manual(values = c(pals[2], pals[1]), name = "Intervention") +
  theme(legend.position = c(0.4, 0.7), text = element_text(size = 12))+
  ylab("Slide prevalence (%) in under 5s")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 60)+
  ggtitle("Additive model")


antag_prev <- ggplot(antag_mod, aes(x = t, y = slide_prev0to5*100, colour = as.factor(int)))+
  geom_line()+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  scale_colour_manual(values = c(pals[2], pals[1]), name = "Intervention") +
  theme(legend.position  = "none", text = element_text(size = 12))+
  ylab("Slide prevalence (%) in under 5s")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 60)+
  ggtitle("Antagonistic model")

plot_prev <- plot_grid(antag_prev, add_prev, ncol = 2,
                       labels = c("C)", "D)"),
                       label_x = 0.035)

dynamics_plot <- plot_grid(EIR_plot, plot_prev, nrow = 2)
ggsave(dynamics_plot, file = "plots/esa-figures/fig_1_dynamics.svg")


add_mod %>%
  group_by(int) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5))
#LLIN only - 40.3
#LLIN + IVM - 8.16


antag_mod %>%
  group_by(int) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5))



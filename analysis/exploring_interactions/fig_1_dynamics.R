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

pals <- c('#1b9e77','#d95f02','#7570b3','#e7298a')

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
#ggsave(dynamics_plot, file = "plots/esa-figures/fig_1_dynamics.svg")


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

antag_mod <- antag_mod %>%
  mutate(model = "antagonistic")

add_mod <- add_mod %>%
  mutate(model = "additive")

both_models <- rbind(antag_mod, add_mod) %>%
  mutate(model_int =  case_when(model == "antagonistic" & int == "LLIN" ~ "antag_LLIN",
                                model == "antagonistic" & int == "LLIN & IVM" ~ "antag_LLIN_IVM",
                                model == "additive" &     int == "LLIN" ~ "add_LLIN",
                                model == "additive" &     int == "LLIN & IVM" ~ "add_LLIN_IVM",
                                TRUE ~ NA_character_))

eir_plot <- ggplot(both_models, aes(x = t, y = EIR_tot, colour = as.factor(model_int)))+
  geom_line()+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  #facet_wrap(~int)+
  scale_colour_manual(values = pals, name = "Model & Intervention",
                      labels = c("Additive: LLIN only", "Additive: LLIN & endectocide",
                                 "Antagonistic: LLIN only", "Antagonistic: LLIN & endectocide")) +
  theme(legend.position  = c(0.3, 0.7), text = element_text(size = 12))+
  ylab("Annual EIR")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 60, yend = 55, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 60, yend = 55, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 60, yend = 55, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 60, yend = 55, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 60, yend = 55, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 75)

prev_plot <- ggplot(both_models, aes(x = t, y = slide_prev0to5*100, colour = as.factor(model_int)))+
  geom_line()+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  #facet_wrap(~int)+
  scale_colour_manual(values = pals, name = "Model & Intervention",
                      labels = c("Additive: LLIN only", "Additive: LLIN & endectocide",
                                 "Antagonistic: LLIN only", "Antagonisitc: LLIN & endectocide")) +
  theme(legend.position  = "none", text = element_text(size = 12))+
  ylab("Slide prevalence (%) in under 5-year-olds")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 45, yend = 41, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 45, yend = 41, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 45, yend = 41, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 45, yend = 41, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 45, yend = 41, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 50)

dynamics_plot2 <- plot_grid(eir_plot, prev_plot, nrow = 2)
ggsave(dynamics_plot2, file = "plots/esa-figures/fig_1_dynamics_plot2.svg")

#modifying plot to put antagonistic and additive model on the same but facet by intervention

#LLIN and IVM
#LLIN only
#IVM only

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

antag_mod_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_antag_IVM.csv", header = TRUE) %>%
  filter(itn_cov == 0 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov, EIR_tot, slide_prev0to5) %>%
  mutate(model = "antagonistic", int = "IVM", id = row_number())

add_mod_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_add_IVM.csv", header = TRUE) %>%
  filter(itn_cov == 0 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov,EIR_tot, slide_prev0to5) %>%
  mutate(model = "additive", int = "IVM", id = row_number())

antag_mod_no_int <- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE) %>%
  filter(itn_cov == 0 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov, EIR_tot, slide_prev0to5) %>%
  mutate(model = "antagonistic", int = "no int", id = row_number())

add_mod_no_int <- read.csv("analysis/esa-analysis/output/pyr_out_add.csv", header = TRUE) %>%
  filter(itn_cov == 0 & d_ITN0 ==  res[2], !t %in% rm) %>%
  select(t, mv, d_ITN0, itn_cov, EIR_tot, slide_prev0to5) %>%
  mutate(model = "antagonistic", int = "no int", id = row_number())

add_mod <- do.call("rbind", list(add_mod_LLIN, add_mod_LLIN_IVM, add_mod_IVM, add_mod_no_int))
antag_mod <- do.call("rbind", list(antag_mod_LLIN, antag_mod_LLIN_IVM, antag_mod_IVM, antag_mod_no_int))

#then rbind the two models
models <- rbind(add_mod, antag_mod)

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
pals <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
          "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")



prev_plot_int_models <- ggplot(models, aes(x = t, y = slide_prev0to5*100))+
  geom_line(aes(colour = as.factor(model)))+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  scale_colour_manual(values = c(pals[2], pals[1]), name = "model") +
  theme(legend.position  = "none", text = element_text(size = 12))+
  ylab("Slide prevalence (%) in under 5s")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 100)+
  facet_wrap(~int)

eir_plot_int_models <- ggplot(models, aes(x = t, y = EIR_tot, colour = as.factor(model)))+
  geom_line()+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)+
  theme_minimal()+
  scale_colour_manual(values = c(pals[2], pals[1]), name = "model") +
  theme(legend.position = c(0.1, 0.2), text = element_text(size = 12))+
  ylab("EIR")+
  scale_x_continuous(limits = c(2250, 3650), breaks = breaks,
                     name = "Time (days)")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 50, yend = 40, colour = "black", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[3], xend = net_seq[3], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  annotate("segment", x = net_seq[4], xend = net_seq[4], y = 55, yend = 40, colour = "blue", arrow = arrow(length = unit(0.02, "npc")))+
  ylim(0, 150)+
  facet_wrap(~int)

plot_int_models <- plot_grid(prev_plot_int_models, eir_plot_int_models,
                             nrow = 2)

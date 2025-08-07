require(tidyverse)

#models runs in net_age_runs_poster.R

antag_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/antag.rds")
add_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/add.rds")

mods <- rbind(antag_mod, add_mod)

res <- unique(antag_mod$d_ITN0)
itn_cov <- unique(antag_mod$itn_cov)

mods_dynamics <- mods %>%
  filter(d_ITN0 == res[1] & model %in% c("antag_LLIN", "antag_LLIN_IVM") & itn_cov == 0.8)

itn_on <- 100 #introduce nets 100 days into simulation



net_seq <- seq(100, 3650, by = 3*365)

#ivm on when nets are 6 months old
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution
mda_int <- 30
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

#when nets are 1yo

IVM_begin2 <- net_seq[3]+(1*365)
IVM_start2 <- c(IVM_begin2, IVM_begin2+mda_int, IVM_begin2 + mda_int + mda_int)

y2.5 <- (365*2.5)
#when nets are 2.5yo
IVM_begin3 <- net_seq[3]+y2.5
IVM_start3 <- c(IVM_begin3, IVM_begin3+mda_int, IVM_begin3+mda_int+mda_int)

lines <- c("antag_LLIN" = "dotted", "antag_LLIN_IVM" = "solid")

sp_pals <- c('#1b9e77','#d95f02','#7570b3','#e7298a')

#DYNAMICS

dynamics_plot <- ggplot(mods_dynamics, aes(x = t/365, y = slide_prev0to5*100, linetype = as.factor(model), col = as.factor(species)))+
  geom_line(linewidth = 1)+
  theme_minimal()+
  labs(y = "Slide prevalence in under 5s (%)",
       x = "Time (years)")+
  scale_linetype_manual(values = lines, name = "Intervention", labels = c("LLIN", "LLIN & endectocide"))+
  ylim(0, 60)+
  annotate("segment", x = IVM_start1[1]/365, xend = IVM_start1[1]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[2]/365, xend = IVM_start1[2]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[3]/365, xend = IVM_start1[3]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 45, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[3]/365, xend = net_seq[3]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[4]/365, xend = net_seq[4]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = (net_seq[1]/365)+2, y = 47, label = "LLIN distribution (every 3 years)", col = "#1f78b4", size = 6)+
  annotate(geom = "text", x = (IVM_start1[1]/365)+(420/365), y = 36, label = "Endectocide MDA", size = 6)+
  scale_colour_manual(name = "Anopheles species", labels = c("arabiensis",
                                                             "funestus",
                                                             "gambiae",
                                                             "stephensi"),
                      values = sp_pals)+
  theme(legend.position = c(0.5, 0.9), legend.direction = "horizontal")+
  scale_x_continuous(breaks=seq(0, 10, 2))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  annotate("rect", fill = "yellow", alpha = 0.5,
           xmin = IVM_start1[1]/365, xmax = (IVM_start1[3]/365)+(23/365),
           ymin = -Inf, ymax = 30)

ggsave(dynamics_plot, file = "analysis/exploring_interactions/MIM_poster/dynamics_plot.svg",
       width = 23.7,
       height = 19.57,
       units = "cm")

#proportion of cases averted when adding endectocides

rel_red_prev <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  filter(itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, species, slide_prev0to5, model) %>%
  spread(key = model, value = slide_prev0to5) %>%
  group_by(species, d_ITN0) %>%
  summarise(mean_prev_antag_LLIN = mean(antag_LLIN),
            mean_prev_antag_LLIN_IVM = mean(antag_LLIN_IVM)) %>%
  mutate(red_prev = ((mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM)/mean_prev_antag_LLIN)*100,
         abs_prev = mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM)

res_pals <- c("#1f78b4", "#33a02c")

no_res <- rel_red_prev %>%
  filter(d_ITN0 == res[1])



efficacy_species_plot <- ggplot(no_res, aes(x = species, y = red_prev))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Relative reduction in slide prevalence in under 5s \n due to endectocide")+
  xlab("Anopheles species")

ggsave(efficacy_species_plot, file = "analysis/exploring_interactions/MIM_poster/eff_sp_plot.svg")

#absolute drop is very very small
efficacy_species_plot_abs <- ggplot(no_res, aes(x = species, y = abs_prev))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Absolute reduction in slide prevalence in under 5s \n due to endectocide")+
  xlab("Anopheles species")

ggsave(efficacy_species_plot_abs, file = "analysis/exploring_interactions/MIM_poster/eff_sp_plot_abs.svg")


#by EIR
rel_red_eir <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  filter(itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, species, EIR_tot, model) %>%
  spread(key = model, value = EIR_tot) %>%
  group_by(species, d_ITN0) %>%
  summarise(mean_eir_antag_LLIN = mean(antag_LLIN),
            mean_eir_antag_LLIN_IVM = mean(antag_LLIN_IVM)) %>%
  mutate(red_eir = ((mean_eir_antag_LLIN - mean_eir_antag_LLIN_IVM)/mean_eir_antag_LLIN)*100,
         abs_eir = mean_eir_antag_LLIN - mean_eir_antag_LLIN_IVM)

res_pals <- c("#1f78b4", "#33a02c")

no_res_eir <- rel_red_eir %>%
  filter(d_ITN0 == res[1])



efficacy_species_plot_abs_eir <- ggplot(no_res_eir, aes(x = species, y = abs_eir))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Absolute reduction in EIR due to endectocide")+
  xlab("Anopheles species")

mod_dynamics_res <- mods %>%
  filter(model %in% c("antag_LLIN", "antag_LLIN_IVM") & itn_cov == 0.8)

dynamics_res_plot <- ggplot(mod_dynamics_res, aes(x = t/365, y = slide_prev0to5*100, linetype = as.factor(model), col = as.factor(species)))+
  geom_line(linewidth = 1)+
  theme_minimal()+
  labs(y = "Slide prevalence in under 5s (%)",
       x = "Time (years)")+
  scale_linetype_manual(values = lines, name = "Intervention", labels = c("LLIN", "LLIN & endectocide"))+
  ylim(0, 48)+
  annotate("segment", x = IVM_start1[1]/365, xend = IVM_start1[1]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[2]/365, xend = IVM_start1[2]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[3]/365, xend = IVM_start1[3]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 45, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[3]/365, xend = net_seq[3]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[4]/365, xend = net_seq[4]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = (net_seq[1]/365)+1, y = 47, label = "LLIN distribution (every 3 years)", col = "#1f78b4")+
  annotate(geom = "text", x = (IVM_start1[1]/365)+(420/365), y = 36, label = "Endectocide MDA campaign")+
  scale_colour_manual(name = "Anopheles species", labels = c("arabiensis",
                                                             "funestus",
                                                             "gambiae",
                                                             "stephensi"),
                      values = sp_pals)+
  theme(legend.position = c(0.35, 0.15), legend.direction = "horizontal")+
  facet_wrap(vars(d_ITN0))

#model performance plot

model_compare_prev <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  #filter(itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, itn_cov, species, slide_prev0to5, model) %>%
  spread(key = model, value = slide_prev0to5) %>%
  group_by(species, d_ITN0, itn_cov) %>%
  summarise(mean_prev_antag_LLIN = mean(antag_LLIN),
            mean_prev_antag_LLIN_IVM = mean(antag_LLIN_IVM),
            mean_prev_add_LLIN = mean(add_LLIN),
            mean_prev_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_prev_antag = ((mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM)/mean_prev_antag_LLIN)*100,
         red_prev_add = ((mean_prev_add_LLIN - mean_prev_add_LLIN_IVM)/mean_prev_add_LLIN)*100,
         abs_prev_antag = mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM,
         abs_prev_add = mean_prev_add_LLIN - mean_prev_add_LLIN_IVM)


model_compare_EIR <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  #filter(itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, itn_cov, species, EIR_tot, model) %>%
  spread(key = model, value = EIR_tot) %>%
  group_by(species, d_ITN0, itn_cov) %>%
  summarise(mean_EIR_antag_LLIN = mean(antag_LLIN),
            mean_EIR_antag_LLIN_IVM = mean(antag_LLIN_IVM),
            mean_EIR_add_LLIN = mean(add_LLIN),
            mean_EIR_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_EIR_antag = ((mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM)/mean_EIR_antag_LLIN)*100,
         red_EIR_add = ((mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)/mean_EIR_add_LLIN)*100,
         abs_EIR_antag = mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM,
         abs_EIR_add = mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)


lm_models <- lm(red_prev_add ~ red_prev_antag, data = model_compare_prev)
summary(lm_models) #adj r2 0.9659
adjusted_r_squared <- summary(lm_models)$adj.r.squared

model_comparison_plot_prev <- ggplot(model_compare_prev, aes(x = red_prev_antag, y = red_prev_add, col = species, shape = as.factor(d_ITN0), size = as.factor(itn_cov),
                          group = 1))+
  geom_point(alpha = 0.8)+
  scale_colour_manual(values = sp_pals, labels = c("arabiensis", "funestus", "gambiae", "stephensi"),
                      name = "Anopheles species")+
  scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
                     name = "Pyrethroid resistance")+
  labs(x = "Efficacy of endectocide (model A)", y = "Efficacy of endectocide (model B)")+
  theme_minimal()+
  theme(legend.position = c(0.9, 0.3))+
  ylim(3, 8)+
  xlim(3, 8)+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed")+
  annotate("text", label = "italic(R)^2 == 0.9659", parse = TRUE, x =6, y = 7, col = "blue", size = 6)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

ggsave(model_comparison_plot, file = "analysis/exploring_interactions/MIM_poster/model_compare.svg",
       width = 23.71,
       height = 16.8,
       units = "cm")
#
lm_models <- lm(abs_EIR_add ~ abs_EIR_antag, data = model_compare_EIR)
summary(lm_models) #adj r2 0.9962
adjusted_r_squared <- summary(lm_models)$adj.r.squared

model_compare_EIR_plot <- ggplot(model_compare_EIR, aes(x = abs_EIR_antag, y = abs_EIR_add, col = species, shape = as.factor(d_ITN0), size = as.factor(itn_cov),
                                                             group = 1))+
  geom_point(alpha = 0.8)+
  scale_colour_manual(values = sp_pals, labels = c("arabiensis", "funestus", "gambiae", "stephensi"),
                      name = "Anopheles species")+
  scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
                     name = "Pyrethroid resistance")+
  labs(x = "Abs red EIR due to endectocide (model A)", y = "Abs red EIR due to endectocide (model B)")+
  theme_minimal()+
  theme(legend.position = c(0.8, 0.3))+
  ylim(0, 6)+
  xlim(0, 6)+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed")+
  annotate("text", label = "italic(R)^2 == 0.9962", parse = TRUE, x =2, y = 4, col = "blue", size = 6)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))




##
ggplot(model_compare, aes(x = red_prev_antag, y = red_prev_add))+
  geom_point()

antag_SM <- readRDS("analysis/exploring_interactions/MIM_poster/antag_SM.rds")
add_SM <- readRDS("analysis/exploring_interactions/MIM_poster/add_SM.rds")
mods_SM <- rbind(antag_SM, add_SM)

model_compare_SM <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  #filter(itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, itn_cov, species, slide_prev0to5, model) %>%
  spread(key = model, value = slide_prev0to5) %>%
  group_by(species, d_ITN0, itn_cov) %>%
  summarise(mean_prev_antag_LLIN = mean(antag_LLIN),
            mean_prev_antag_LLIN_IVM = mean(antag_LLIN_IVM),
            mean_prev_add_LLIN = mean(add_LLIN),
            mean_prev_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_prev_antag = ((mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM)/mean_prev_antag_LLIN)*100,
         red_prev_add = ((mean_prev_add_LLIN - mean_prev_add_LLIN_IVM)/mean_prev_add_LLIN)*100)

ggplot(model_compare_SM, aes(x = red_prev_antag, y = red_prev_add, col = species, shape = as.factor(d_ITN0), size = as.factor(itn_cov)))+
  geom_point()


ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")

ivm_haz <- as.data.frame(ivm_haz)

hr_plot <- ggplot(ivm_haz, aes(x = Day, y = IVM_300_3_HS))+
  geom_point(size = 2)+
  labs(x = expression("Day of bloodmeal after first dose of ivermectin-like drug (3x300 " * mu * "g/kg)"),
       y = "Hazard ratio")+
  theme_minimal()+
  geom_vline(xintercept = 23, lty = "dashed", col = "black")+
  ylim(0, 10)+
  theme(text = element_text(size = 22))+
  geom_hline(yintercept = 1, lty = "dashed", col = "red")

ggsave(hr_plot, file = "analysis/exploring_interactions/MIM_poster/hazards.svg",
       width = 30.31,
       height = 14.22,
       units = "cm")


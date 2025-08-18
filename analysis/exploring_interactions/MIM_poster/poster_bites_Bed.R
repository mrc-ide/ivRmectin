require(tidyverse)

#model runs in bites_Bed analysis

antag_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/bites_Bed/antag.rds")
add_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/bites_Bed/add.rds")

mods <- rbind(antag_mod, add_mod)
#NB this is 10% resistance.
res <- unique(antag_mod$d_ITN0)
itn_cov <- unique(antag_mod$itn_cov)

mods_dynamics <- mods %>%
  mutate(model_type = case_when(model %in% c("antag_LLIN", "antag_LLIN_IVM") ~ "antag",
                                 model %in% c("add_LLIN", "add_LLIN_IVM") ~ "add",
                                 TRUE ~ NA_character_))

itn_on <- 100 #introduce nets 100 days into simulation



net_seq <- seq(100, 3650, by = 3*365)
mda_int <- 30
#ivm on when nets are 6 months old
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution
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
bb_labs <- unique(mods_dynamics$bites_Bed)
dynamics_plot <- ggplot(mods_dynamics, aes(x = t/365, y = slide_prev0to5*100, linetype = as.factor(model), col = as.factor(bites_Bed)))+
  geom_line(linewidth = 1)+
  theme_minimal()+
  labs(y = "Slide prevalence in under 5s (%)",
       x = "Time (years)")+
  scale_linetype_manual(values = lines, name = "Intervention", labels = c("LLIN", "LLIN & endectocide"))+
  ylim(0, 100)+
  annotate("segment", x = IVM_start1[1]/365, xend = IVM_start1[1]/365, y = 65, yend = 61, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[2]/365, xend = IVM_start1[2]/365, y = 65, yend = 61, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[3]/365, xend = IVM_start1[3]/365, y = 65, yend = 61, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 70, yend = 65, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 65, yend = 61, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[3]/365, xend = net_seq[3]/365, y = 65, yend = 61, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[4]/365, xend = net_seq[4]/365, y = 65, yend = 61, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = (net_seq[1]/365)+3, y = 75, label = "LLIN distribution (every 3 years)", col = "#1f78b4", size = 6)+
  annotate(geom = "text", x = (IVM_start1[1]/365)+(420/365), y = 68, label = "Endectocide MDA", size = 6)+
  scale_colour_manual(name = "Proportion of bites in bed", labels = bb_labs,
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
           ymin = -Inf, ymax = 61)

ggsave(dynamics_plot, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/dynamics_plot.svg",
       width = 23.7,
       height = 19.57,
       units = "cm")

#proportion of cases averted when adding endectocides
##ANTAG MODELS
rel_red_prev <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  filter(model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, bites_Bed, slide_prev0to5, model) %>%
  spread(key = model, value = slide_prev0to5) %>%
  group_by(bites_Bed) %>%
  summarise(mean_prev_antag_LLIN = mean(antag_LLIN),
            mean_prev_antag_LLIN_IVM = mean(antag_LLIN_IVM)) %>%
  mutate(red_prev = ((mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM)/mean_prev_antag_LLIN)*100,
         abs_prev = mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM)

res_pals <- c("#1f78b4", "#33a02c")


efficacy_species_plot_rel_prev <- ggplot(rel_red_prev, aes(x = as.factor(bites_Bed), y = red_prev))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Relative reduction in slide prevalence in under 5s \n due to endectocide")+
  xlab("Proportion of bites in bed")
ggsave(efficacy_species_plot_rel_prev, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/eff_sp_plot_rel_prev.svg")

efficacy_species_plot_abs_prev <- ggplot(rel_red_prev, aes(x = as.factor(bites_Bed), y = abs_prev))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Absolute reduction in slide prevalence in under 5s \n due to endectocide")+
  xlab("Proportion of bites in bed")+
  ylim(0, 0.03)

ggsave(efficacy_species_plot_abs_prev, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/eff_sp_plot_abs_prev.svg")

rel_red_EIR <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  filter(model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, bites_Bed, EIR_tot, model) %>%
  spread(key = model, value = EIR_tot) %>%
  group_by(bites_Bed) %>%
  summarise(mean_EIR_antag_LLIN = mean(antag_LLIN),
            mean_EIR_antag_LLIN_IVM = mean(antag_LLIN_IVM)) %>%
  mutate(red_EIR = ((mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM)/mean_EIR_antag_LLIN)*100,
         abs_EIR = mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM)

efficacy_species_plot_rel_EIR <- ggplot(rel_red_EIR, aes(x = as.factor(bites_Bed), y = red_EIR))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Relative reduction in EIR \n due to endectocide")+
  xlab("Proportion of bites in bed")
ggsave(efficacy_species_plot_rel_EIR, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/eff_sp_plot_rel_EIR.svg")

efficacy_species_plot_abs_EIR <- ggplot(rel_red_EIR, aes(x = as.factor(bites_Bed), y = abs_EIR))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Absolute reduction in EIR in under 5s \n due to endectocide")+
  xlab("Proportion of bites in bed")

ggsave(efficacy_species_plot_abs_EIR, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/eff_sp_plot_abs_EIR.svg")

#efficacy for the additive models
red_prev_add <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  filter(model %in% c("add_LLIN", "add_LLIN_IVM")) %>%
  select(t, d_ITN0, bites_Bed, slide_prev0to5, model) %>%
  spread(key = model, value = slide_prev0to5) %>%
  group_by(bites_Bed) %>%
  summarise(mean_prev_add_LLIN = mean(add_LLIN),
            mean_prev_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_prev = ((mean_prev_add_LLIN - mean_prev_add_LLIN_IVM)/mean_prev_add_LLIN)*100,
         abs_prev = mean_prev_add_LLIN - mean_prev_add_LLIN_IVM)

red_EIR_add <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  filter(model %in% c("add_LLIN", "add_LLIN_IVM")) %>%
  select(t, d_ITN0, bites_Bed, EIR_tot, model) %>%
  spread(key = model, value = EIR_tot) %>%
  group_by(bites_Bed) %>%
  summarise(mean_EIR_add_LLIN = mean(add_LLIN),
            mean_EIR_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_EIR = ((mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)/mean_EIR_add_LLIN)*100,
         abs_EIR = mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)

#look at absolute drops in EIR and prevalence with additive model
#with high phi-B, nets kill before mozzie has had the chance to take up ivermectin, meaning endectocides are less effective

ggplot(red_prev_add, aes(x = as.factor(bites_Bed), y = abs_prev))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Absolute reduction in prev in under 5s \n due to endectocide")+
  xlab("Proportion of bites in bed")

ggplot(red_EIR_add, aes(x = as.factor(bites_Bed), y = abs_EIR))+
  geom_bar(stat = "identity", fill = "white", col = "black")+
  theme_minimal()+
  #scale_colour_manual(values = res_pals, labels = c("90% resistance", "No resistance"),
  #                    name = "Pyrethroid resistance")+
  ylab("Absolute reduction in EIR \n due to endectocide")+
  xlab("Proportion of bites in bed")



##
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

model_compare <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  #filter(itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, itn_cov, bites_Bed, slide_prev0to5, model) %>%
  spread(key = model, value = slide_prev0to5) %>%
  group_by(bites_Bed, d_ITN0, itn_cov) %>%
  summarise(mean_prev_antag_LLIN = mean(antag_LLIN),
            mean_prev_antag_LLIN_IVM = mean(antag_LLIN_IVM),
            mean_prev_add_LLIN = mean(add_LLIN),
            mean_prev_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(rel_red_prev_antag = ((mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM)/mean_prev_antag_LLIN)*100,
         rel_red_prev_add = ((mean_prev_add_LLIN - mean_prev_add_LLIN_IVM)/mean_prev_add_LLIN)*100,
         abs_red_prev_antag = mean_prev_antag_LLIN - mean_prev_antag_LLIN_IVM,
         abs_red_prev_add = mean_prev_add_LLIN - mean_prev_add_LLIN_IVM)


lm_models <- lm(rel_red_prev_add ~ rel_red_prev_antag, data = model_compare)
summary(lm_models) #adj r2 0.9598
adjusted_r_squared <- summary(lm_models)$adj.r.squared

model_comparison_plot_rel_prev <- ggplot(model_compare, aes(x = rel_red_prev_antag, y = rel_red_prev_add, col = as.factor(bites_Bed), group = 1))+
  geom_point(alpha = 0.8)+
  scale_colour_manual(values = sp_pals, labels = c("0.25", "0.5", "0.75", "0.9"),
                      name = "Proportion of bites in bed")+
  #scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  #scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
  #                   name = "Pyrethroid resistance")+
  labs(x = "Efficacy of endectocide (model A) antag", y = "Efficacy of endectocide (model B) add")+
  theme_minimal()+
  theme(legend.position = c(0.7, 0.8))+
  ylim(3, 8)+
  xlim(3, 8)+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed")+
  annotate("text", label = "italic(R)^2 == 0.9598", parse = TRUE, x =9, y = 7, col = "blue", size = 6)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

summary(lm(abs_red_prev_add ~ abs_red_prev_antag, data = model_compare))

model_comparison_plot_abs_prev <- ggplot(model_compare, aes(x = abs_red_prev_antag, y = abs_red_prev_add, col = as.factor(bites_Bed), group = 1))+
  geom_point(alpha = 0.8, size = 5)+
  scale_colour_manual(values = sp_pals, labels = c("0.25", "0.5", "0.75", "0.9"),
                      name = "Proportion of bites in bed")+
  #scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  #scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
  #                   name = "Pyrethroid resistance")+
  labs(x = "Abs red prev due to endectocide (model A) antag", y = "Abs red prev due to endectocide (model B) add")+
  theme_minimal()+
  theme(legend.position = c(0.7, 0.8))+
  ylim(0, 0.05)+
  xlim(0, 0.05)+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed")+
  annotate("text", label = "italic(R)^2 == 0.9598", parse = TRUE, x =9, y = 7, col = "blue", size = 6)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))



ggsave(model_comparison_plot, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/model_compare.svg",
       width = 23.71,
       height = 16.8,
       units = "cm")

model_compare_EIR <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  #filter(d_ITN0 == res[1] & itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  #filter(itn_cov == 0.8 & model %in% c("antag_LLIN", "antag_LLIN_IVM")) %>%
  select(t, d_ITN0, itn_cov, bites_Bed, EIR_tot, model) %>%
  spread(key = model, value = EIR_tot) %>%
  group_by(bites_Bed, d_ITN0, itn_cov) %>%
  summarise(mean_EIR_antag_LLIN = mean(antag_LLIN),
            mean_EIR_antag_LLIN_IVM = mean(antag_LLIN_IVM),
            mean_EIR_add_LLIN = mean(add_LLIN),
            mean_EIR_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_EIR_antag = ((mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM)/mean_EIR_antag_LLIN)*100,
         red_EIR_add = ((mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)/mean_EIR_add_LLIN)*100,
         abs_EIR_antag = mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM,
         abs_EIR_add = mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)

summary(lm(abs_EIR_add ~ abs_EIR_antag, data = model_compare_EIR)) #R2 is 0.999

lm_models <- lm(red_EIR_add ~ red_EIR_antag, data = model_compare_EIR)
summary(lm_models) #adj r2 0.9574
adjusted_r_squared <- summary(lm_models)$adj.r.squared

model_comparison_plot_EIR <- ggplot(model_compare_EIR, aes(x = abs_EIR_antag, y = abs_EIR_add, col = as.factor(bites_Bed), group = 1))+
  geom_point(alpha = 0.8)+
  scale_colour_manual(values = sp_pals, labels = c("0.25", "0.5", "0.75", "0.9"),
                      name = "Proportion of bites in bed")+
  #scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  #scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
  #                   name = "Pyrethroid resistance")+
  labs(x = "Absolute reduction in EIR due to endectocide (model A)", y = "Absolute reduction in EIR due to endectocide (model B)")+
  theme_minimal()+
  theme(legend.position = c(0.7, 0.8))+
  ylim(0, 3)+
  xlim(0, 3)+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed", col = "grey")+
  annotate("text", label = "italic(R)^2 == 0.999", parse = TRUE, x =2.5, y = 2, col = "red", size = 6)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave(model_comparison_plot, file = "analysis/exploring_interactions/MIM_poster/bites_Bed/model_compare.svg",
       width = 23.71,
       height = 16.8,
       units = "cm")

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


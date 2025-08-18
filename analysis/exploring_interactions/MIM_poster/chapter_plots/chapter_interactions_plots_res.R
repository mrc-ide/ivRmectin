#plotting bites bed, Q0, res and coverage
require(tidyverse)
require(RColorBrewer)
#dynamics (prevalence) then absolute reduction in by both models


#initial setup
itn_on <- 100 #introduce nets 100 days into simulation

net_seq <- seq(100, 3650, by = 3*365)
mda_int <- 30
#ivm on when nets are 6 months old
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

bb_pal <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
Q0_pal <- c('#fdcc8a','#fc8d59','#e34a33','#b30000')
res_pal <- c('#b3cde3','#8c96c6','#8856a7','#810f7c')
cov_pal <- c('#d7b5d8','#df65b0','#dd1c77','#980043')
sp_pals <- c('#1b9e77','#d95f02','#7570b3','#e7298a')
lines <- c("antag_LLIN" = "dotted", "antag_LLIN_IVM" = "solid")

#res plots####
antag_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/itn_res/antag.rds")
add_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/itn_res/add.rds")
head(add_mod)
mods <- rbind(antag_mod, add_mod)
unique(mods$Q0)
res <- unique(antag_mod$d_ITN0)
itn_cov <- unique(antag_mod$itn_cov)

mods_dynamics <- mods %>%
  mutate(model_type = case_when(model %in% c("antag_LLIN", "antag_LLIN_IVM") ~ "antag",
                                model %in% c("add_LLIN", "add_LLIN_IVM") ~ "add",
                                TRUE ~ NA_character_))


#DYNAMICS
res_labs <- unique(mods_dynamics$d_ITN0)
res_labs2 <- c("80", "60", "20", "0")
res_name <- "Death in mortality bioassay (%)"
dynamics_plot_res <- ggplot(antag_mod, aes(x = t/365, y = slide_prev0to5*100, linetype = as.factor(model), col = as.factor(d_ITN0)))+
  geom_line(linewidth = 1)+
  theme_minimal()+
  labs(y = "Slide prevalence in under 5s (%) (model A)",
       x = "Time (years)")+
  scale_linetype_manual(values = lines, name = "Intervention", labels = c("LLIN", "LLIN & endectocide"))+
  ylim(0, 100)+
  annotate("segment", x = IVM_start1[1]/365, xend = IVM_start1[1]/365, y = 50, yend = 41, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[2]/365, xend = IVM_start1[2]/365, y = 50, yend = 41, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[3]/365, xend = IVM_start1[3]/365, y = 50, yend = 41, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 50, yend = 42, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 50, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[3]/365, xend = net_seq[3]/365, y = 50, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[4]/365, xend = net_seq[4]/365, y = 50, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = (net_seq[1]/365)+3, y = 55, label = "LLIN distribution (every 3 years)", col = "#1f78b4", size = 3)+
  annotate(geom = "text", x = (IVM_start1[1]/365)+(420/365), y = 52, label = "Endectocide MDA", size = 3)+
  scale_colour_manual(name = "d_ITN0", values = res_pal)+
  theme(legend.position = c(0.5, 0.8), legend.direction = "vertical")+
  scale_x_continuous(breaks=seq(0, 10, 2))+
  #theme(axis.text.x = element_text(size = 12),
  #      axis.text.y = element_text(size = 12),
  #      axis.title.x = element_text(size = 12),
  #      axis.title.y = element_text(size = 12),
  #      legend.text = element_text(size = 12),
  #      legend.title = element_text(size = 12))+
  annotate("rect", fill = "yellow", alpha = 0.5,
           xmin = IVM_start1[1]/365, xmax = (IVM_start1[3]/365)+(23/365),
           ymin = -Inf, ymax = 40)

model_compare_EIR <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, d_ITN0, model, EIR_tot, itn_cov) %>%
  spread(key = model, value = EIR_tot) %>%
  group_by(d_ITN0) %>%
  summarise(mean_EIR_antag_LLIN = mean(antag_LLIN),
            mean_EIR_antag_LLIN_IVM = mean(antag_LLIN_IVM),
            mean_EIR_add_LLIN = mean(add_LLIN),
            mean_EIR_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_EIR_antag = ((mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM)/mean_EIR_antag_LLIN)*100,
         red_EIR_add = ((mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)/mean_EIR_add_LLIN)*100,
         abs_EIR_antag = mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM,
         abs_EIR_add = mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)

summary(lm(abs_EIR_add ~ abs_EIR_antag, data = model_compare_EIR)) #R2 is ~1

model_comparison_plot_abs_EIR_res <- ggplot(model_compare_EIR, aes(x = abs_EIR_antag, y = abs_EIR_add, col = as.factor(d_ITN0), group = 1))+
  geom_point(size = 3)+
  scale_color_manual(values = res_pal, labels = res_labs2,
                      name = "Phenotypic insecticide resistance")+
  #scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  #scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
  #                   name = "Pyrethroid resistance")+
  labs(x = "Absolute reduction in EIR due to endectocide (model A)", y = "Absolute reduction in EIR due to endectocide (model B)")+
  theme_minimal()+
  #theme(legend.position = c(0.8, 0.8))+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed", col = "grey", alpha = 0.5)+
  annotate("text", label = "italic(R)^2 == 0.999", parse = TRUE, x = 3, y = 4, col = "red", size = 6)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  guides(color = "none")+
  xlim(0, 5)+
  ylim(0, 5)

res_traits <- antag_mod %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, avhc, mu, d_ITN0) %>%
  group_by(d_ITN0) %>%
  summarise(av_time_bloodmeals = 1/mean(avhc), av_life_exp = 1/mean(mu))

res_traits_antag <- gather(res_traits, trait, value, av_time_bloodmeals:av_life_exp,
                          factor_key = TRUE) %>%
  mutate(model = "antag")

res_traits_add <- add_mod %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, avhc, mu, d_ITN0) %>%
  group_by(d_ITN0) %>%
  summarise(av_time_bloodmeals = 1/mean(avhc), av_life_exp = 1/mean(mu))

res_traits_add <- gather(res_traits_add, trait, value, av_time_bloodmeals:av_life_exp,
                        factor_key = TRUE) %>%
  mutate(model = "add")

res_traits <- rbind(res_traits_antag, res_traits_add)

model_labels <- c(
  `av_time_bloodmeals` = "Average time \n between human bloodmeals",
  `av_life_exp` = "Average mosquito \n life expectancy")

traits_plot <- ggplot(res_traits, aes(x = as.factor(d_ITN0), y = value, col = as.factor(model)))+
  geom_bar(stat = "identity", position = position_dodge(), fill = "white", size = 1)+
  #scale_x_discrete(labels = c("Average time \n between human bloodmeals",
  #                            "Average mosquito \n life expectancy"), name = "Mosquito life trait (model B)")+
  ylab("Time (days)")+
  theme_minimal()+
  #scale_fill_manual(values = bb_pal, labels = c("0.25", "0.5", "0.75", "0.9"),
  #                  name = "Proportion of bites in bed")+
  guides(fill = "none")+
  facet_wrap(vars(trait), labeller = as_labeller(model_labels))+
  scale_colour_manual(labels = c("model B", "model A"), values = c("#d8b365", '#5ab4ac'),
                      name = "Model")+
  ylim(0, 8)+
  xlab("d_ITN0")+
  theme(legend.position = c(0.3, 0.8))+
  scale_x_discrete(labels = c("0.22", "0.27", "0.32", "0.34"))

res_plots <- cowplot::plot_grid(dynamics_plot_res, model_comparison_plot_abs_EIR_res, traits_plot, labels = c("A)", "B)", "C)"), ncol = 3)
ggsave(res_plots, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/res_plots.svg")





#itn_cov####

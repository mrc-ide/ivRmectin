#initial set up####
#initial setup
require(tidyverse)
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

# Q0####
antag_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/Q0/antag.rds")
add_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/Q0/add.rds")

mods <- rbind(antag_mod, add_mod)
Q0_labs <- unique(antag_mod$Q0)

mods_dynamics <- mods %>%
  mutate(model_type = case_when(model %in% c("antag_LLIN", "antag_LLIN_IVM") ~ "antag",
                                model %in% c("add_LLIN", "add_LLIN_IVM") ~ "add",
                                TRUE ~ NA_character_))

#DYNAMICS
dynamics_plot_Q0 <- ggplot(antag_mod, aes(x = t/365, y = slide_prev0to5*100, linetype = as.factor(model), col = as.factor(Q0)))+
  geom_line(linewidth = 1)+
  theme_minimal()+
  labs(y = "Slide prevalence in under 5s (%) (model A)",
       x = "Time (years)")+
  scale_linetype_manual(values = lines, name = "Intervention", labels = c("LLIN", "LLIN & endectocide"))+
  ylim(0, 100)+
  annotate("segment", x = IVM_start1[1]/365, xend = IVM_start1[1]/365, y = 25, yend = 21, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[2]/365, xend = IVM_start1[2]/365, y = 25, yend = 21, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[3]/365, xend = IVM_start1[3]/365, y = 25, yend = 21, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 45, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 45, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[3]/365, xend = net_seq[3]/365, y = 45, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[4]/365, xend = net_seq[4]/365, y = 45, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = (net_seq[1]/365)+3, y = 46, label = "LLIN distribution (every 3 years)", col = "#1f78b4", size = 3)+
  annotate(geom = "text", x = (IVM_start1[1]/365)+(420/365), y = 26, label = "Endectocide MDA", size = 3)+
  scale_colour_manual(name = "Human blood index", labels = Q0_labs, values = Q0_pal)+
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
           ymin = -Inf, ymax = 20)

model_compare_EIR <- mods %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, d_ITN0, itn_cov, Q0, EIR_tot, model) %>%
  spread(key = model, value = EIR_tot) %>%
  group_by(Q0, d_ITN0, itn_cov) %>%
  summarise(mean_EIR_antag_LLIN = mean(antag_LLIN),
            mean_EIR_antag_LLIN_IVM = mean(antag_LLIN_IVM),
            mean_EIR_add_LLIN = mean(add_LLIN),
            mean_EIR_add_LLIN_IVM = mean(add_LLIN_IVM)) %>%
  mutate(red_EIR_antag = ((mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM)/mean_EIR_antag_LLIN)*100,
         red_EIR_add = ((mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)/mean_EIR_add_LLIN)*100,
         abs_EIR_antag = mean_EIR_antag_LLIN - mean_EIR_antag_LLIN_IVM,
         abs_EIR_add = mean_EIR_add_LLIN - mean_EIR_add_LLIN_IVM)


summary(lm(abs_EIR_add ~ abs_EIR_antag, data = model_compare_EIR)) #R2 is 0.981

model_comparison_plot_abs_EIR_Q0 <- ggplot(model_compare_EIR, aes(x = abs_EIR_antag, y = abs_EIR_add, col = as.factor(Q0), group = 1))+
  geom_point(size = 3)+
  scale_colour_manual(values = Q0_pal, labels = c("0.25", "0.5", "0.75", "0.9"),
                      name = "Human Blood Index")+
  #scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  #scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
  #                   name = "Pyrethroid resistance")+
  labs(x = "Absolute reduction in EIR due to endectocide (model A)", y = "Absolute reduction in EIR due to endectocide (model B)")+
  theme_minimal()+
  #theme(legend.position = c(0.8, 0.8))+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed", col = "grey", alpha = 0.5)+
  annotate("text", label = "italic(R)^2 == 0.981", parse = TRUE, x = 2, y = 3, col = "red", size = 6)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  guides(color = "none")+
  xlim(0, 3)+
  ylim(0, 3)

bb_traits <- antag_mod %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, avhc, mu, Q0) %>%
  group_by(Q0) %>%
  summarise(av_time_bloodmeals = 1/mean(avhc), av_life_exp = 1/mean(mu))

bb_traits_antag <- gather(bb_traits, trait, value, av_time_bloodmeals:av_life_exp,
                          factor_key = TRUE) %>%
  mutate(model = "antag")

bb_traits_add <- add_mod %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, avhc, mu, Q0) %>%
  group_by(Q0) %>%
  summarise(av_time_bloodmeals = 1/mean(avhc), av_life_exp = 1/mean(mu))

bb_traits_add <- gather(bb_traits_add, trait, value, av_time_bloodmeals:av_life_exp,
                        factor_key = TRUE) %>%
  mutate(model = "add")

bb_traits <- rbind(bb_traits_antag, bb_traits_add)

model_labels <- c(
  `av_time_bloodmeals` = "Average time \n between human bloodmeals",
  `av_life_exp` = "Average mosquito \n life expectancy")

traits_plot_Q0 <- ggplot(bb_traits, aes(x = Q0, y = value, col = as.factor(model)))+
  geom_bar(stat = "identity", position = "dodge", fill = "white", size = 1)+
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
  xlab("Human Blood Index")+
  ylim(0, 32)+
  theme(legend.position = c(0.3, 0.8))

Q0_plots <- cowplot::plot_grid(dynamics_plot_Q0, model_comparison_plot_abs_EIR_Q0, traits_plot_Q0,
                   labels  = c("A)", "B)", "C)"), ncol = 3)

ggsave(Q0_plots, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/Q0_plots.svg")

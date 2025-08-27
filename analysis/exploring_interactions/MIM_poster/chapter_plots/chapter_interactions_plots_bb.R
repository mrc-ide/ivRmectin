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

#bites bed plots####
antag_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/bites_Bed/antag.rds")
add_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/bites_Bed/add.rds")
head(add_mod)
mods <- rbind(antag_mod, add_mod)
res <- unique(antag_mod$d_ITN0)
itn_cov <- unique(antag_mod$itn_cov)

mods_dynamics <- mods %>%
  mutate(model_type = case_when(model %in% c("antag_LLIN", "antag_LLIN_IVM") ~ "antag",
                                model %in% c("add_LLIN", "add_LLIN_IVM") ~ "add",
                                TRUE ~ NA_character_))



mods_dynamics %>%
  filter(model_type == "add") %>%
  ggplot(aes(x = t/365, y = slide_prev0to5*100))+
  geom_line()

#DYNAMICS
bb_labs <- unique(mods_dynamics$bites_Bed)
dynamics_plot_bb <- ggplot(antag_mod, aes(x = t/365, y = slide_prev0to5*100, linetype = as.factor(model), col = as.factor(bites_Bed)))+
  geom_line(linewidth = 1)+
  theme_minimal()+
  labs(y = "Slide prevalence in under 5-year-olds (%) \n (model A)",
       x = "Time (years)")+
  scale_linetype_manual(values = lines, name = "Intervention", labels = c("ITN", "ITN & endectocide"))+
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
  scale_colour_manual(name = "Proportion of bites \n when people are in bed", labels = bb_labs, values = bb_pal)+
  theme(legend.position = c(0.5, 0.8), legend.direction = "vertical")+
  scale_x_continuous(breaks=seq(0, 10, 2))+
  guides(color = "none")
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
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  mutate(EIRout = EIRout*1000) %>% #per 1000 persons
  select(t, d_ITN0, itn_cov, bites_Bed, EIRout, model) %>%
  spread(key = model, value = EIRout) %>%
  group_by(bites_Bed, d_ITN0, itn_cov) %>%
  summarise(tot_EIR_antag_LLIN = sum(antag_LLIN),
            tot_EIR_antag_LLIN_IVM = sum(antag_LLIN_IVM),
            tot_EIR_add_LLIN = sum(add_LLIN),
            tot_EIR_add_LLIN_IVM = sum(add_LLIN_IVM)) %>%
  mutate(red_EIR_antag = ((tot_EIR_antag_LLIN - tot_EIR_antag_LLIN_IVM)/tot_EIR_antag_LLIN)*100,
         red_EIR_add = ((tot_EIR_add_LLIN - tot_EIR_add_LLIN_IVM)/tot_EIR_add_LLIN)*100,
         abs_EIR_antag = tot_EIR_antag_LLIN - tot_EIR_antag_LLIN_IVM,
         abs_EIR_add = tot_EIR_add_LLIN - tot_EIR_add_LLIN_IVM)

EIR_stats <- summary(lm(abs_EIR_add ~ abs_EIR_antag, data = model_compare_EIR)) #R2 is 0.9997
EIR_Rsq <- round(EIR_stats$adj.r.squared, 3)

model_comparison_plot_abs_EIR_bb <- ggplot(model_compare_EIR, aes(x = abs_EIR_antag, y = abs_EIR_add, fill = as.factor(bites_Bed), group = 1))+
  geom_point(size = 3, shape = 21, colour = "black", stroke = 1) +  # black outline
  scale_fill_manual(values = bb_pal, labels = c("0.25", "0.5", "0.75", "0.9"),
                      name = "Proportion of bites taken on people \n when they are in bed")+
  #scale_size_manual(values = c(3, 5), labels = c("20%", "80%"), name = "LLIN coverage")+
  #scale_shape_manual(values = c(19, 17), labels = c("90% resistance", "No resistance"),
  #                   name = "Pyrethroid resistance")+
  labs(x = "Absolute reduction in EIR per 1000 persons \n due to endectocide (model A)", y = "Absolute reduction in EIR per 1000 persons \n due to endectocide (model B)")+
  theme_minimal()+
  theme(legend.position = c(0.4, 0.6))+
  geom_smooth(method="lm", se = FALSE, show.legend = FALSE, lty = "dashed", col = "grey", alpha = 0.5)+
  annotate("text",
           label = paste0("italic(R)^2 == ", EIR_Rsq),
           parse = TRUE,
           x = 20,
           y = 2500,
           col = "red", size = 6, hjust = 0)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) #+
  #guides(color = "none")+
  #coord_cartesian(xlim = c(0, 15), ylim = c(0, 15))

bb_traits <- antag_mod %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, avhc, mu, bites_Bed) %>%
  group_by(bites_Bed) %>%
  summarise(av_time_bloodmeals = 1/mean(avhc), av_life_exp = 1/mean(mu))

bb_traits_antag <- gather(bb_traits, trait, value, av_time_bloodmeals:av_life_exp,
                         factor_key = TRUE) %>%
  mutate(model = "antag")

bb_traits_add <- add_mod %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  select(t, avhc, mu, bites_Bed) %>%
  group_by(bites_Bed) %>%
  summarise(av_time_bloodmeals = 1/mean(avhc), av_life_exp = 1/mean(mu))

bb_traits_add <- gather(bb_traits_add, trait, value, av_time_bloodmeals:av_life_exp,
                         factor_key = TRUE) %>%
  mutate(model = "add")

bb_traits <- rbind(bb_traits_antag, bb_traits_add)

model_labels <- c(
  `av_time_bloodmeals` = "Average time \n between human bloodmeals",
  `av_life_exp` = "Average mosquito \n life expectancy")

traits_plot <- ggplot(bb_traits, aes(x = bites_Bed, y = value, fill = as.factor(model)))+
  geom_bar(stat = "identity", position = "dodge", col = "black", size = 1)+
  #scale_x_discrete(labels = c("Average time \n between human bloodmeals",
  #                            "Average mosquito \n life expectancy"), name = "Mosquito life trait (model B)")+
  ylab("Time (days)")+
  theme_minimal()+
  #scale_fill_manual(values = bb_pal, labels = c("0.25", "0.5", "0.75", "0.9"),
  #                  name = "Proportion of bites in bed")+
  facet_wrap(vars(trait), labeller = as_labeller(model_labels))+
  scale_fill_manual(labels = c("model B", "model A"), values = c("#d8b365", '#5ab4ac'),
                      name = "Model")+
  ylim(0, 8)+
  xlab("Proportion of bites \n when people are in bed")+
  theme(legend.position = c(0.3, 0.8))

bb_plots <- cowplot::plot_grid(dynamics_plot_bb, model_comparison_plot_abs_EIR_bb, traits_plot, labels = c("A)", "B)", "C)"), ncol = 3)

#ggsave(bb_plots, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/bb_plots.svg")


#then do the heatmaps for each combination of bites_Bed and other parameters
#in bites_Bed_heatpmaps.R (copied here)
require(tidyverse)

df_bb_Q0_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN_IVM.rds")
df_bb_cov_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN_IVM.rds")
df_bb_res_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN_IVM.rds")

df_bb_Q0_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN.rds")
df_bb_cov_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN.rds")
df_bb_res_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN.rds")

#lower Q0 files
df_bb_cov_ITN_IVM2 <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_2_ITN_IVM.rds")
df_bb_res_ITN_IVM2 <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_2_ITN_IVM.rds")

df_bb_cov_ITN2 <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_2_ITN.rds")
df_bb_res_ITN2 <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_2_ITN.rds")


df_bb_cov_ITN_plot <- df_bb_cov_ITN %>%
  mutate(IVRM_sr =0)

df_bb_cov_plot <- rbind(df_bb_cov_ITN_plot,df_bb_cov_ITN_IVM)

lines_int <- c("ITN" = "dotted", "ITN_IVM" = "solid")

dynamics_plot_bb_odin <- df_bb_cov_plot %>%
  filter(itn_cov == 0.8) %>%
  ggplot()+
  aes(x = t/365, y = slide_prev0to5*100, lty = as.factor(int), col = as.factor(bites_Bed))+
  geom_line(size = 1)+
  theme_bw()+
  scale_linetype_manual(values = lines_int, name = "Intervention", labels = c("ITN", "ITN & endectocide"))+
  ylim(0, 100)+
  annotate("segment", x = IVM_start1[1]/365, xend = IVM_start1[1]/365, y = 70, yend = 60, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[2]/365, xend = IVM_start1[2]/365, y = 70, yend = 60, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[3]/365, xend = IVM_start1[3]/365, y = 70, yend = 60, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 75, yend = 65, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 70, yend = 60, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[3]/365, xend = net_seq[3]/365, y = 70, yend = 60, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[4]/365, xend = net_seq[4]/365, y = 70, yend = 60, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = (net_seq[1]/365)+3, y = 55, label = "ITN distribution (every 3 years)", col = "#1f78b4", size = 3)+
  annotate(geom = "text", x = (IVM_start1[1]/365)+(420/365), y = 52, label = "Endectocide MDA", size = 3)+
  scale_colour_manual(name = "Proportion of bites \n when people are in bed", labels = bb_labs, values = bb_pal)+
  theme(legend.position = c(0.5, 0.8), legend.direction = "vertical")+
  scale_x_continuous(breaks=seq(0, 10, 2))+
  guides(color = "none")


df_bb_cov_plot %>%
  filter(itn_cov == 0.8) %>%
  ggplot()+
  aes(x = t/365, y = clin_inc0to5*1000, lty = as.factor(int), col = as.factor(bites_Bed))+
  geom_line(size = 1)+
  theme_bw()+
  coord_cartesian(xlim = c(6.4, 7.5))+
  scale_linetype_manual(values = lines_int, name = "Intervention", labels = c("ITN", "ITN & endectocide"))+
  geom_vline(xintercept = 6.76)+
  geom_vline(xintercept = 7.3)

df_bb_cov_plot %>%
  filter(itn_cov == 0.8) %>%
  ggplot()+
  aes(x = t/365, y = EIRout*1000, lty = as.factor(int), col = as.factor(bites_Bed))+
  geom_line(size = 1)+
  theme_bw()+
  coord_cartesian(xlim = c(6.4, 7.5))+
  scale_linetype_manual(values = lines_int, name = "Intervention", labels = c("ITN", "ITN & endectocide"))+
  geom_vline(xintercept = 6.76)+
  geom_vline(xintercept = 7.3)

df_bb_cov_plot %>%
  filter(itn_cov == 0.8) %>%
  ggplot()+
  aes(x = t/365, y = mv, lty = as.factor(int), col = as.factor(bites_Bed))+
  geom_line(size = 1)+
  theme_bw()+
  coord_cartesian(xlim = c(6.4, 7.5))+
  scale_linetype_manual(values = lines_int, name = "Intervention", labels = c("ITN", "ITN & endectocide"))+
  geom_vline(xintercept = 6.76)+
  geom_vline(xintercept = 7.3)



net_seq <- seq(100, 3650, by = 3*365)
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution

#common limits for plots
common_limits_rel <- c(20, 45)


#heatmap bites_Bed and Q0

output_bb_Q0_ITN_IVM <- df_bb_Q0_ITN_IVM %>%
  select(t, bites_Bed, Q0, clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  rename(clin_inc_bb_Q0_ITN_IVM = clin_inc0to5)


output_bb_Q0_ITN <- df_bb_Q0_ITN %>%
  select(t, bites_Bed, Q0, clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  rename(clin_inc_ITN = clin_inc0to5)

output_bb_Q0 <- left_join(output_bb_Q0_ITN, output_bb_Q0_ITN_IVM)

output_bb_Q0 <- output_bb_Q0 %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(Q0, bites_Bed) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_bb_Q0_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100,
         species = case_when(Q0 == 0.95 & bites_Bed == 0.95 ~ "italic('An. gambiae')",
                             Q0 == 0.95 & bites_Bed == 0.75 ~ "italic('An. funestus')",
                             Q0 == 0.75 & bites_Bed == 0.75 ~ "italic('An. arabiensis')",
                             Q0 == 0.25 & bites_Bed == 0.50 ~ "italic('An. stephensi')",
                             TRUE ~ ""))

heatmap_Q0_bb_rel <- ggplot(output_bb_Q0, aes(x = as.factor(bites_Bed), y = as.factor(Q0), fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits_rel, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Human Blood Index")+
  guides(fill = "none")+
  #geom_text(aes(label = species), parse = TRUE, col = "white", size = 5)+
  geom_text(aes(label = paste0(round(rel_diff_IVM, 1), "%")),
            col = "white", size = 5)


common_limits <- c(50, 350)

heatmap_Q0_bb_abs <- ggplot(output_bb_Q0, aes(x = as.factor(bites_Bed), y = as.factor(Q0), fill = abs_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Human Blood Index")+
  guides(fill = "none")+
  #geom_text(aes(label = species), parse = TRUE, col = "white", size = 5)+
  geom_text(aes(label = paste0(round(abs_diff_IVM))),
            col = "white", size = 5)

ggplot(output_bb_Q0, aes(x = as.factor(bites_Bed), y = as.factor(Q0), fill = tot_cases_ITN_IVM))+
  geom_tile()+
  theme_bw()+
  #scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Human Blood Index")+
  guides(fill = "none")+
  #geom_text(aes(label = species), parse = TRUE, col = "white", size = 5)+
  geom_text(aes(label = paste0(round(tot_cases_ITN_IVM))),
            col = "white", size = 5)+
  ggtitle("ITN & IVM")


ggplot(output_bb_Q0, aes(x = as.factor(bites_Bed), y = as.factor(Q0), fill = tot_cases_ITN))+
  geom_tile()+
  theme_bw()+
  #scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Human Blood Index")+
  guides(fill = "none")+
  #geom_text(aes(label = species), parse = TRUE, col = "white", size = 5)+
  geom_text(aes(label = paste0(round(tot_cases_ITN))),
            col = "white", size = 5)+
  ggtitle("ITN")





#heatmap bites_Bed and cov

output_bb_cov_ITN_IVM <- df_bb_cov_ITN_IVM %>%
  select(t, bites_Bed, itn_cov, clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>% #make it per 1000 persons
  rename(clin_inc_bb_cov_ITN_IVM = clin_inc0to5)


output_bb_cov_ITN <- df_bb_cov_ITN %>%
  select(t, bites_Bed, itn_cov, clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>% #make it per 1000 persons
  rename(clin_inc_ITN = clin_inc0to5)

output_bb_cov <- left_join(output_bb_cov_ITN, output_bb_cov_ITN_IVM)

output_bb_cov <- output_bb_cov %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(bites_Bed, itn_cov) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_bb_cov_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100)

heatmap_bb_cov_rel <- ggplot(output_bb_cov, aes(x = as.factor(bites_Bed), y = as.factor(itn_cov*100), fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits_rel, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("ITN coverage %")+
  guides(fill = "none")+
  geom_text(aes(label = paste0(round(rel_diff_IVM, 1), "%")),col = "white", size = 5)

heatmap_bb_cov_abs <- ggplot(output_bb_cov, aes(x = as.factor(bites_Bed), y = as.factor(itn_cov*100), fill = abs_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("ITN coverage %")+
  guides(fill = "none")+
  geom_text(aes(label = paste0(round(abs_diff_IVM))),col = "white", size = 5)

EIR_ITN_IVM <- df_bb_cov_ITN_IVM %>%
  #filter(itn_cov == 0.8) %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  group_by(itn_cov, bites_Bed) %>%
  summarise(tot_EIR_bb_cov_ITN_IVM = sum(clin_inc0to5))

EIR_ITN <- df_bb_cov_ITN %>%
  #filter(itn_cov == 0.8) %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  group_by(itn_cov, bites_Bed) %>%
  summarise(tot_EIR_bb_cov_ITN= sum(clin_inc0to5))

out_EIR_bb_cov <- left_join(EIR_ITN_IVM, EIR_ITN) %>%
  mutate(rel_diff_IVM = ((tot_EIR_bb_cov_ITN - tot_EIR_bb_cov_ITN_IVM)/tot_EIR_bb_cov_ITN)*100,
         abs_diff_IVM = tot_EIR_bb_cov_ITN - tot_EIR_bb_cov_ITN_IVM)

ggplot(out_EIR_bb_cov, aes(x = as.factor(bites_Bed), y = as.factor(itn_cov*100), fill = abs_diff_IVM))+
  geom_tile()+
  theme_bw()+
  #scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("ITN coverage %")+
  geom_text(aes(label = round(abs_diff_IVM, 2)),col = "white", size = 5)


ggplot(out_EIR_bb_cov, aes(x = as.factor(bites_Bed), y = as.factor(itn_cov*100), fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()+
  #scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("ITN coverage %")+
  geom_text(aes(label = round(rel_diff_IVM, 2)),col = "white", size = 5)




##bb cov2
ggplot(df_bb_cov_ITN_IVM2, aes(x = t, y = EIRout, col = as.factor(itn_cov*100)))+
  geom_line()+
  facet_wrap(vars(bites_Bed))

ggplot(df_bb_cov_ITN_IVM2, aes(x = t, y = slide_prev0to5, col = as.factor(itn_cov*100)))+
  geom_line()+
  facet_wrap(vars(bites_Bed))

ggplot(df_bb_cov_ITN_IVM2, aes(x = t, y = clin_inc0to5, col = as.factor(itn_cov*100)))+
  geom_line()+
  facet_wrap(vars(bites_Bed))

ggplot(df_bb_cov_ITN_IVM2, aes(x = t, y = mv, col = as.factor(itn_cov*100)))+
  geom_line()+
  facet_wrap(vars(bites_Bed))

output_bb_cov2_ITN_IVM <- df_bb_cov_ITN_IVM2 %>%
  select(t, bites_Bed, itn_cov, EIRout) %>%
  rename(EIRout_bb_cov2_ITN_IVM = EIRout)

output_bb_cov2_ITN <- df_bb_cov_ITN2 %>%
  select(t, bites_Bed, itn_cov, EIRout) %>%
  rename(EIRout_bb_cov2_ITN = EIRout)

output_bb_cov2 <- left_join(output_bb_cov2_ITN, output_bb_cov2_ITN_IVM)

check <- output_bb_cov2 %>%
  filter(between(t, IVM_begin1, IVM_begin1+180))

output_bb_cov2 <- output_bb_cov2 %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(bites_Bed, itn_cov) %>%
  summarise(tot_EIR_ITN_IVM = sum(EIRout_bb_cov2_ITN_IVM),
            tot_EIR_ITN = sum(EIRout_bb_cov2_ITN))%>%
  mutate(abs_diff_IVM = tot_EIR_ITN - tot_EIR_ITN_IVM,
         rel_diff_IVM = ((tot_EIR_ITN - tot_EIR_ITN_IVM)/tot_EIR_ITN)*100)

heatmap_bb_cov2 <- ggplot(output_bb_cov2, aes(x = as.factor(bites_Bed), y = as.factor(itn_cov*100), fill = abs_diff_IVM))+
  geom_tile()+
  theme_bw()+
  #scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("ITN coverage %")+
  geom_text(aes(label = round(rel_diff_IVM, 2)),col = "white", size = 5)





#heatmap bites_bed and res
output_bb_res_ITN_IVM <- df_bb_res_ITN_IVM %>%
  select(t, bites_Bed, d_ITN0, r_ITN0,  clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  rename(clin_inc_bb_res_ITN_IVM = clin_inc0to5)


output_bb_res_ITN <- df_bb_res_ITN %>%
  select(t, bites_Bed, d_ITN0, r_ITN0, clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  rename(clin_inc_ITN = clin_inc0to5)

output_bb_res <- left_join(output_bb_res_ITN, output_bb_res_ITN_IVM)

output_bb_res <- output_bb_res %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(bites_Bed, d_ITN0, r_ITN0) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_bb_res_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100)

res_bb_df <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/chapter_interactions_bb_res_df.rds")

res_bb_df1 <- res_bb_df %>%
  select(d_ITN0, r_ITN0, resistance) %>%
  distinct()

output_bb_res <- left_join(output_bb_res, res_bb_df1)


heatmap_bb_res_rel <- ggplot(output_bb_res, aes(x = as.factor(bites_Bed), y = as.factor(resistance*100), fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits_rel, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Phenotypic resistance (%)")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")+
  geom_text(aes(label = paste0(round(rel_diff_IVM, 1), "%")),col = "white", size = 5)

heatmap_bb_res_abs <- ggplot(output_bb_res, aes(x = as.factor(bites_Bed), y = as.factor(resistance*100), fill = abs_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits, name = "Cases averted (absolute) in \n under 5-year-olds due to endectocide",
                       breaks = seq(50, 350, length.out = 5))+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Phenotypic resistance (%)")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.5, "cm"))+
  geom_text(aes(label = round(abs_diff_IVM)),col = "white", size = 5)


fig_1_odin_model_B <- cowplot::plot_grid(dynamics_plot_bb_odin, model_comparison_plot_abs_EIR_bb, traits_plot,
                                       heatmap_Q0_bb_abs, heatmap_bb_cov_abs, heatmap_bb_res_abs,
                                       labels = c("A", "B", "C", "D", "E", "F"), ncol = 3)


fig_1_odin_model_C <- cowplot::plot_grid(dynamics_plot_bb_odin, model_comparison_plot_abs_EIR_bb, traits_plot,
                                         heatmap_Q0_bb_rel, heatmap_bb_cov_rel, heatmap_bb_res_rel,
                                         labels = c("A", "B", "C", "D", "E", "F"), ncol = 3)




fig_1_odin_model <- cowplot::plot_grid(dynamics_plot_bb, model_comparison_plot_abs_EIR_bb, traits_plot,
                                       heatmap_Q0_bb, heatmap_bb_cov, heatmap_bb_res,
                                       labels = c("A", "B", "C", "D", "E", "F"), ncol = 3)

heatmap_Q0_bb_abs

## bb-res2
##bb cov2
output_bb_res2_ITN_IVM <- df_bb_res_ITN_IVM2 %>%
  select(t, bites_Bed, d_ITN0, clin_inc0to5) %>%
  rename(clin_inc_bb_res2_ITN_IVM = clin_inc0to5)

output_bb_res2_ITN <- df_bb_res_ITN2 %>%
  select(t, bites_Bed, d_ITN0, clin_inc0to5) %>%
  rename(clin_inc_bb_res2_ITN = clin_inc0to5)

output_bb_res2 <- left_join(output_bb_res2_ITN, output_bb_res2_ITN_IVM)

output_bb_res2 <- output_bb_res2 %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(bites_Bed, d_ITN0) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_bb_res2_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_bb_res2_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100)

heatmap_bb_res2 <- ggplot(output_bb_res2, aes(x = as.factor(bites_Bed), y = as.factor(d_ITN0), fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("d_ITN0")+
  geom_text(aes(label = round(rel_diff_IVM, 2)),col = "white", size = 5)



##with high bites_Bed, there is little variability in impact due to d_ITN0. Check inputs

ggplot(output_bb_res, aes(x = bites_Bed, y = r_ITN0, fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(limits = common_limits, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Human Blood Index")

fig_1_odin_model <- cowplot::plot_grid(dynamics_plot_bb, model_comparison_plot_abs_EIR_bb, traits_plot,
                  heatmap_Q0_bb, heatmap_bb_cov, heatmap_bb_res,
                  labels = c("A", "B", "C", "D", "E", "F"), ncol = 3)




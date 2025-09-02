#dynamics plots, efficacy plots and heatmaps of influence of different bionomics parameters
require(tidyverse)
df_bb_Q0_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN_IVM.rds")
df_bb_cov_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN_IVM.rds")
df_bb_res_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN_IVM.rds")

df_bb_Q0_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN.rds")
df_bb_cov_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN.rds")
df_bb_res_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN.rds")


#initial setup
itn_on <- 100 #introduce nets 100 days into simulation

net_seq <- seq(100, 3650, by = 3*365)
mda_int <- 30
#ivm on when nets are 6 months old
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)


bb_pal <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

df_bb_cov_ITN_plot <- df_bb_cov_ITN %>%
  mutate(IVRM_sr =0)

df_bb_cov_plot <- rbind(df_bb_cov_ITN_plot,df_bb_cov_ITN_IVM)

lines_int <- c("ITN" = "dotted", "ITN_IVM" = "solid")
bb_labs <- unique(df_bb_cov_plot$bites_Bed)

#net_distr_3 <- (t-net_seq[3])/365
itn_covs <- unique(df_bb_cov_ITN_plot$itn_cov)
dynamics_plot_bb_odin <- df_bb_cov_plot %>%
  filter(itn_cov == itn_covs[3]) %>%
  ggplot()+
  aes(x = (t-net_seq[3])/365, y = slide_prev0to5*100, lty = as.factor(int), col = as.factor(bites_Bed))+
  geom_line(size = 1)+
  theme_bw(base_size = 14)+
  scale_linetype_manual(values = lines_int, name = "Intervention", labels = c("ITN", "ITN & endectocide"))+
  coord_cartesian(xlim = c(0,2.8), ylim = c(0,100))+
  annotate("segment",
           x = (IVM_start1[1] - net_seq[3])/365,
           xend = (IVM_start1[1] - net_seq[3])/365,
           y = 70, yend = 60,
           colour = "black",
           arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment",
        x = (IVM_start1[2] - net_seq[3])/365,
        xend = (IVM_start1[2] - net_seq[3])/365,
        y = 70, yend = 60,
        colour = "black",
        arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment",
           x = (IVM_start1[3] - net_seq[3])/365,
           xend = (IVM_start1[3] - net_seq[3])/365,
           y = 70, yend = 60,
           colour = "black",
           arrow = arrow(length = unit(0.01, "npc")))+
  #annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 75, yend = 65, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  #annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 70, yend = 60, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = 0, xend = 0, y = 70, yend = 60, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = 3, xend = 3, y = 70, yend = 60, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = 1, y = 85, label = "ITN distribution (every 3 years)", col = "#1f78b4")+
  annotate(geom = "text", x = 1.2, y = 75, label = "Endectocide MDA")+
  scale_colour_manual(name = "Proportion of bites \n when people are in bed", labels = bb_labs, values = bb_pal)+
  theme(legend.position = c(0.7, 0.15), legend.direction = "vertical")+
  scale_x_continuous(breaks=seq(0, 10, 2))+
  guides(color = "none")+
  labs(y = "Slide prevalence (%) \n in under 5-year-olds", x = "Time (years) since most recent ITN campaign")

#read in the constant uptake models
cons_uptake_bb_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_IVM_ITN.rds")
cons_uptake_bb_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cons_ITN.rds")

EIR_ITN_IVM_cons <- cons_uptake_bb_ITN_IVM %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180) & Q0 == 0.95) %>%
  mutate(EIRout_1000_persons = EIRout*1000) %>%
  group_by(bites_Bed) %>%
  summarise(EIR_total_ITN_IVM = sum(EIRout_1000_persons))

EIR_ITN_cons <- cons_uptake_bb_ITN %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180) & Q0 == 0.95) %>%
  mutate(EIRout_1000_persons = EIRout*1000) %>%
  group_by(bites_Bed) %>%
  summarise(EIR_total_ITN = sum(EIRout_1000_persons))

EIR_summary_cons <- left_join(EIR_ITN_cons, EIR_ITN_IVM_cons) %>%
  mutate(abs_diff_EIR_cons = (EIR_total_ITN - EIR_total_ITN_IVM),
         rel_diff_EIR_cons = ((EIR_total_ITN - EIR_total_ITN_IVM)/EIR_total_ITN)*100)

#antag summary output

antag_uptake_bb_ITN_IVM <- df_bb_cov_ITN_IVM %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  filter(itn_cov == itn_covs[3]) %>%
  mutate(EIRout_1000_persons = EIRout*1000) %>%
  group_by(bites_Bed) %>%
  summarise(EIR_total_ITN_IVM = sum(EIRout_1000_persons))

antag_uptake_bb_ITN <- df_bb_cov_ITN %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  filter(itn_cov == itn_covs[3]) %>%
  mutate(EIRout_1000_persons = EIRout*1000) %>%
  group_by(bites_Bed) %>%
  summarise(EIR_total_ITN= sum(EIRout_1000_persons))

EIR_summary_antag <- left_join(antag_uptake_bb_ITN,antag_uptake_bb_ITN_IVM)%>%
  mutate(abs_diff_EIR_antag = (EIR_total_ITN - EIR_total_ITN_IVM),
         rel_diff_EIR_antag = ((EIR_total_ITN - EIR_total_ITN_IVM)/EIR_total_ITN)*100)

EIR_summary_cons <- EIR_summary_cons %>%
  select(bites_Bed, abs_diff_EIR_cons, rel_diff_EIR_cons)

EIR_summary_antag <- EIR_summary_antag %>%
  select(bites_Bed, abs_diff_EIR_antag, rel_diff_EIR_antag)

mod_compare_EIR <- left_join(EIR_summary_antag, EIR_summary_cons)

#pearsons corr coef
EIR_stats <- with(mod_compare_EIR, cor(rel_diff_EIR_antag, rel_diff_EIR_cons)) #correlation strength


EIR_efficacy_plot <- ggplot(mod_compare_EIR, aes(x = rel_diff_EIR_antag, y = rel_diff_EIR_cons, fill = as.factor(bites_Bed)))+
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 50))+
  #xlim(0,50)+
  #ylim(0,50)
  theme_bw(base_size = 14)+
  geom_point(size = 3, shape = 21, colour = "black", stroke = 1) +  # black outline
  scale_fill_manual(values = bb_pal, labels = c("0.25", "0.5", "0.75", "0.9"),
                    name = expression(phi[italic(Bed)]))+
  theme(legend.position = c(0.7, 0.2))+
  labs(x = "Model A: Efficacy in reducing EIR (%)",y = "Model B: Efficacy in reducing EIR (%)")+
  annotate("text",
           label = paste0("italic(r) == ", round(EIR_stats,2)),
           parse = TRUE,
           x = 15,
           y = 40,
           col = "blue", size = 6, hjust = 0)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40")  # y = x line

#avhc and mu plots
#need to bring in the traits from the constant uptake model

Q0_vals <- unique(df_bb_Q0_ITN$Q0)
unique(df_bb_Q0_ITN$itn_cov)
traits <- df_bb_Q0_ITN %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  #filter(itn_cov == 0.6 & Q0 == 0.95) %>%
  filter(Q0 != Q0_vals[2]) %>%
  group_by(bites_Bed, Q0) %>%
  summarise(mean_time_between_meals = 1/mean(avhc),
            mean_life_exp = 1/mean(mu))

traits_summary <- traits %>%
  pivot_wider(
    id_cols = bites_Bed,
    names_from = Q0,
    values_from = c(mean_time_between_meals, mean_life_exp),
    names_glue = "{.value}_{Q0}"
  ) %>%
  rename(
    mean_inv_avhc = mean_time_between_meals_0.75,
    mean_inv_avhc_lower = mean_time_between_meals_0.25,
    mean_inv_avhc_upper = mean_time_between_meals_0.95,
    mean_inv_mort = mean_life_exp_0.75,
    mean_inv_mort_lower = mean_life_exp_0.25,
    mean_inv_mort_upper = mean_life_exp_0.95
  ) %>%
  mutate(model = "antag")


traits_cons <- cons_uptake_bb_ITN %>%
  filter(between(t, IVM_start1[1], IVM_start1[3]+23)) %>%
  filter(Q0 != Q0_vals[2]) %>%
  group_by(bites_Bed, Q0) %>%
  summarise(mean_time_between_meals = 1/mean(avhc),
            mean_life_exp = 1/mean(mu))

traits_summary_cons <- traits_cons %>%
  pivot_wider(
    id_cols = bites_Bed,
    names_from = Q0,
    values_from = c(mean_time_between_meals, mean_life_exp),
    names_glue = "{.value}_{Q0}"
  ) %>%
  rename(
    mean_inv_avhc = mean_time_between_meals_0.75,
    mean_inv_avhc_lower = mean_time_between_meals_0.25,
    mean_inv_avhc_upper = mean_time_between_meals_0.95,
    mean_inv_mort = mean_life_exp_0.75,
    mean_inv_mort_lower = mean_life_exp_0.25,
    mean_inv_mort_upper = mean_life_exp_0.95
  ) %>%
  mutate(model = "constant-uptake")


traits_combo <- rbind(traits_summary, traits_summary_cons)

model_pals <- c('#7fc97f','#beaed4')

traits_avhc <- ggplot(traits_combo, aes(x = as.factor(bites_Bed), y = mean_inv_avhc, fill = as.factor(model)))+
  geom_bar(stat = "identity", position = position_dodge(0.9), col = "black")+
  geom_errorbar(aes(ymin = mean_inv_avhc_lower, ymax = mean_inv_avhc_upper), position = position_dodge(0.9),
                width = 0.4, size = 1)+
  theme_bw(base_size = 14)+
  labs(x = "Proportion of bites in bed", y = "Average time between \n human bloodmeals (days)")+
  theme(legend.position = c(0.5, 0.85))+
  scale_fill_manual(values = model_pals, labels = c("ITN-mediated endectocide uptake (model A)",
                                                    "Constant endectocide uptake (model B)"),
                    name = "Model assumption")+
  coord_cartesian(ylim = c(0,30))

traits_mort <- ggplot(traits_combo, aes(x = as.factor(bites_Bed), y = mean_inv_mort, fill = as.factor(model)))+
  geom_bar(stat = "identity", position = position_dodge(0.9), col = "black")+
  geom_errorbar(aes(ymin = mean_inv_mort_lower, ymax = mean_inv_mort_upper),
                position = position_dodge(0.9),
                width = 0.4, size = 1)+
  theme_bw(base_size = 14)+
  labs(x = "Proportion of bites in bed", y = "Average mosquito \n life expectancy (days)")+
  guides(fill = "none")+
  scale_fill_manual(values = model_pals, labels = c("ITN-mediated endectocide uptake",
                                                    "Constant endectocide uptake"))


#heatmaps####

df_bb_Q0_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN_IVM.rds")
df_bb_cov_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN_IVM.rds")
df_bb_res_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN_IVM.rds")

df_bb_Q0_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN.rds")
df_bb_cov_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN.rds")
df_bb_res_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN.rds")

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
  theme_bw(base_size = 14)+
  scale_fill_viridis_c(limits = common_limits_rel, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Proportion of bites taken on people when they are in bed")+
  ylab("Human Blood Index")+
  #guides(fill = "none")+
  #geom_text(aes(label = species), parse = TRUE, col = "white", size = 5)+
  geom_text(aes(label = paste0(round(rel_diff_IVM, 1), "%")),
            col = "white", size = 5)+
  theme(legend.position = "bottom")


figure_1_odin_1 <- cowplot::plot_grid(dynamics_plot_bb_odin, EIR_efficacy_plot,
                                    traits_avhc, traits_mort,
                                    labels = c("A", "B", "C", "D"))



figure_1_odin <- cowplot::plot_grid(figure_1_odin_1, heatmap_Q0_bb_rel,
                                    labels = c("", "E"))

ggsave(figure_1_odin, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/figure_1_odin.pdf")
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


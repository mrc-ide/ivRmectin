
#viz and epi impact of staggered vs all in one distribution

#sensitivity analysis: endemicity, (init_EIR = 2 or 100), Q0 (four levels, for each species) and seasonality (no seasonality or Fattick, Senegal)
require(tidyverse)
require(svglite)

start <- (365*5)+200

df_distr_all <- readRDS("analysis/target-profiles-distrib-strat/chapter/output/df_distr_HR_Sen_3m_HS.rds")

model_types <-  unique(df_distr_all$model_type)

df_distr_all %>%
  filter(t == 1) %>% #EIR at t = 1 either 1.39 or 69.5
  group_by(ref) %>%
  summarise(eir = EIR_tot) %>%
  distinct() %>%  #so ref 1-3 are low EIR and ref >3 are high EIR
  print(n = 14)

df_distr_all <- df_distr_all %>%
  mutate(init_EIR = case_when(ref == 1 ~ 2,
                              TRUE ~ 100))

unique(df_distr_all$init_EIR)
unique(df_distr_all$model_type)

df_distr_all %>%
  filter(model_type == "baseline" & init_EIR == 100 & Q0 == 0.92) %>%
  ggplot()+
  aes(x = t, y = clin_inc0to5*1000)+
  geom_line()+
  ylim(0, 6)
#blue is baseline
#red - all in one
#green - 10 days
#purple - 30 days

distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')

#filter to just high EIR and gamb-like vector (Q0 = 0.92) for plots

covs <- unique(df_distr_all$ivm_cov)

df_distr <- df_distr_all %>%
  filter(init_EIR == 100)


#measurements taken at diff time points for trials
#MATAMAL: 4 weeks after last MDA (prevalence)
#BOHEMIA: incidence from first MDA, for 6 months
matamal_survey <- (30*4) + start
plot_matamal <- matamal_survey-start #diff between start and survey
2145-start
bohemia_inc_period <- start+(6*30)
plot_bohemia <- bohemia_inc_period- start

#check pop has been dividided up correctly, can see stagger

#stag_10 <- 10
#times_10_stag <- seq(0, 0+ stag_10, length.out = num_distr)
#
#
#stag_30 <- 30
#times_30_stag <- seq(0, 0+ stag_30, length.out = num_distr)

#by 0.28y since int, the mosquito density is back to equilibrium
mv_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = mv, col = as.factor(model_type))) +
  geom_line(size = 1.1) +
  theme_bw() +
  ylim(0, 300) +
  ylab("Mosquito density")+
  theme(legend.position = c(0.4, 0.8)) +

  scale_color_manual(name = "Scenario", values = distr_pals,
                     labels = c("10 days to complete 1 round of MDA", "20 days to complete 1 round of MDA",
                                "1 day to complete 1 round of MDA", "Baseline")) +
  xlim(-0.25, 1) +
  xlab("Years since intervention started") +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 1, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = plot_matamal/365, col = "grey", linetype = "dashed", size = 1.05)  +
  geom_vline(xintercept = plot_bohemia/365, col = "orange", linetype = "dashed", size = 1.05)

#  # Arrows added with annotate()
#  #all-in-one
#  annotate("segment", x = 0, y = 47, xend = 0, yend = 44,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[1], size = 1) +
#
#  #10 day staggered
#  annotate("segment", x = times_10_stag[2]/365, y = 47, xend = times_10_stag[2]/365, yend = 44,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  annotate("segment", x = times_10_stag[3]/365, y = 47, xend = times_10_stag[3]/365, yend = 44,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  #30d staggered
#  annotate("segment", x = times_30_stag[2]/365, y = 47, xend = times_30_stag[2]/365, yend = 44,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1) +
#
#  annotate("segment", x = times_30_stag[3]/365, y = 47, xend = times_30_stag[3]/365, yend = 44,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1)


#v_count <- df_distr %>%
# filter(ivm_cov == covs[2]) %>%
# select(t, mv, model_type) %>%
# pivot_wider(values_from = mv, names_from = model_type)

#v_count_baseline <- df_distr_all %>%
# filter(model_type == "baseline" & init_EIR == 100 & Q0 == 0.92) %>%
# select(t, mv) %>%
# rename(baseline_mv = mv)

#v_killed_int_stats <- left_join(mv_count, mv_count_baseline) %>%
# mutate(killed_10_d_stagger = baseline_mv - `10d-stagger`,
#        killed_20d_d_stagger = baseline_mv - `20d-stagger`,
#        killed_all_in_one = baseline_mv -  `all-in-one`,
#        tot_10d = cumsum(killed_10_d_stagger),
#        tot_20d = cumsum(killed_20d_d_stagger),
#        tot_all = cumsum(killed_20d_d_stagger))

#umsum(mv_killed_int_stats$killed_10_d_stagger)
#umsum(mv_killed_int_stats$killed_20d_d_stagger)
#umsum(mv_killed_int_stats$killed_all_in_one)

#um(mv_killed_int_stats$killed_10_d_stagger)
#um(mv_killed_int_stats$killed_20d_d_stagger)
#um(mv_killed_int_stats$killed_all_in_one)


#v_killed_int <- mv_killed_int_stats %>%
# select(t, killed_10_d_stagger, killed_20d_d_stagger, killed_all_in_one) %>%
# pivot_longer(cols = -t,names_to = "model_type", values_to = "mv_killed")



#osq_killed <- ggplot(mv_killed_int, aes(x = (t - start)/365, y = mv_killed, col = as.factor(model_type))) +
# geom_line(size = 1.1) +
# theme_bw() +
# theme(legend.position = c(0.8, 0.2)) +

# scale_color_manual(name = "Scenario", values = distr_pals) +
# xlim(-0.25, 1)+
# guides(col = "none")+
# ylab("Number of mosquitoes killed by intervention")+
# xlab("Years since intervention started")+
# geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
# geom_vline(xintercept = 1, col = "black", linetype = "dashed", size = 1.05) +
# geom_vline(xintercept = plot_matamal/365, col = "grey", linetype = "dashed", size = 1.05)  +
# geom_vline(xintercept = plot_bohemia/365, col = "orange", linetype = "dashed", size = 1.05)

# Arrows added with annotate()
#all-in-one
#  annotate("segment", x = 0, y = 30, xend = 0, yend = 26,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[1], size = 1) +
#
#  #10 day staggered
#  annotate("segment", x = times_10_stag[2]/365, y = 30, xend = times_10_stag[2]/365, yend = 26,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  annotate("segment", x = times_10_stag[3]/365, y = 30, xend = times_10_stag[3]/365, yend = 26,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  #30d staggered
#  annotate("segment", x = times_30_stag[2]/365, y = 30, xend = times_30_stag[2]/365, yend = 26,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1) +
#
#  annotate("segment", x = times_30_stag[3]/365, y = 30, xend = times_30_stag[3]/365, yend = 26,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1)

eir_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = EIR_tot, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  ylim(0, 300)+
  scale_color_manual(name = "Scenario", values = distr_pals)+
  guides(col = "none")+
  xlim(-0.25, 1) +
  ylab("Average annual EIR")+
  xlab("Years since intervention started") +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 1, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = plot_matamal/365, col = "grey", linetype = "dashed", size = 1.05)  +
  geom_vline(xintercept = plot_bohemia/365, col = "orange", linetype = "dashed", size = 1.05)

# Arrows added with annotate()
#all-in-one
#  annotate("segment", x = 0, y = 81, xend = 0, yend = 77,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[1], size = 1) +
#
#  #10 day staggered
#  annotate("segment", x = times_10_stag[2]/365, y = 81, xend = times_10_stag[2]/365, yend = 77,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  annotate("segment", x = times_10_stag[3]/365, y = 81, xend = times_10_stag[3]/365, yend = 77,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  #30d staggered
#  annotate("segment", x = times_30_stag[2]/365, y = 81, xend = times_30_stag[2]/365, yend = 77,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1) +
#
#  annotate("segment", x = times_30_stag[3]/365, y = 81, xend = times_30_stag[3]/365, yend = 77,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1)

prev_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = slide_prev0to5*100, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  ylim(0, 80)+
  scale_color_manual(name = "Scenario", values = distr_pals)+
  guides(col = "none")+
  xlim(-0.25, 1) +
  xlab("Years since intervention started")+
  ylab("Slide prevalence (%) in children \n under 5-years-old")+
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 1, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = plot_matamal/365, col = "grey", linetype = "dashed", size = 1.05)  +
  geom_vline(xintercept = plot_bohemia/365, col = "orange", linetype = "dashed", size = 1.05)

#  # Arrows added with annotate()
#  #all-in-one
#  annotate("segment", x = 0, y = 69, xend = 0, yend = 65,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[1], size = 1) +
#
#  #10 day staggered
#  annotate("segment", x = times_10_stag[2]/365, y = 69, xend = times_10_stag[2]/365, yend = 65,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  annotate("segment", x = times_10_stag[3]/365, y = 69, xend = times_10_stag[3]/365, yend = 65,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  #30d staggered
#  annotate("segment", x = times_30_stag[2]/365, y = 69, xend = times_30_stag[2]/365, yend = 65,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1) +
#
#  annotate("segment", x = times_30_stag[3]/365, y = 69, xend = times_30_stag[3]/365, yend = 65,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1)

inc_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = clin_inc0to5*1000, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  scale_color_manual(name = "Scenario", values = distr_pals)+
  guides(col = "none")+
  ylab("Clinical incidence in children \n under 5-years-old, per 1000 persons")+
  ylim(0, 25)+
  xlim(-0.25, 1) +
  xlab("Years since intervention started")+
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 1, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = plot_matamal/365, col = "grey", linetype = "dashed", size = 1.05)  +
  geom_vline(xintercept = plot_bohemia/365, col = "orange", linetype = "dashed", size = 1.05)

#  # Arrows added with annotate()
#  #all-in-one
#  annotate("segment", x = 0, y = 5.4, xend = 0, yend = 5,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[1], size = 1) +
#
#  #10 day staggered
#  annotate("segment", x = times_10_stag[2]/365, y = 5.4, xend = times_10_stag[2]/365, yend = 5,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  annotate("segment", x = times_10_stag[3]/365, y = 5.4, xend = times_10_stag[3]/365, yend = 5,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[3], size = 1) +
#
#  #30d staggered
#  annotate("segment", x = times_30_stag[2]/365, y = 5.4, xend = times_30_stag[2]/365, yend = 5,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1) +
#
#  annotate("segment", x = times_30_stag[3]/365, y = 5.4, xend = times_30_stag[3]/365, yend = 5,
#           arrow = arrow(length = unit(0.3, "cm")),
#           color = distr_pals[4], size = 1)

dynamics <- cowplot::plot_grid(mv_plot, eir_plot, inc_plot, prev_plot, align = "v",
                               labels = c("A", "B", "C", "D"))

#then look at epi impact
#df_int <- do.call("rbind", list(df_mod1, df_mod2, df_mod3))

#by day 300, all back to eqm for mv

#look at epi impact for all

#need to look at this over different time periods. Note, that this is over 5y so overall impact is very diluted.
#slower strategy takes longer to return to eqm hence appears more impactful

#make df_distr_all according to each scenario

df_distr_wide_setting_ivm <- df_distr_all %>%
  select(t,ref, init_EIR, model_type, clin_inc0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)
  }) %>%
  select(-`baseline-Sen`)

df_distr_wide_setting_baseline <- df_distr_all %>%
  filter(model_type == "baseline-Sen") %>%
  select(t,  ref,init_EIR, model_type, clin_inc0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)
  })


df_distr_wide_setting <- left_join(df_distr_wide_setting_ivm,df_distr_wide_setting_baseline, by = c("ref",
                                                                                                    "t",

                                                                                                    "init_EIR"))
df_distr_wide_setting_prev_ivm<- df_distr_all %>%
  select(t, ref, init_EIR, model_type, slide_prev0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to5)
  }) %>%
  select(-`baseline-Sen`)

df_distr_wide_setting_prev_baseline <- df_distr_all %>%
  filter(model_type == "baseline-Sen") %>%
  select(t,ref,init_EIR, model_type, slide_prev0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to5)
  })

df_distr_wide_setting_prev <- left_join(df_distr_wide_setting_prev_ivm,df_distr_wide_setting_prev_baseline, by = c("ref",
                                                                                                                   "t",

                                                                                                                   "init_EIR"))


#measure impact from start to 1y later
#timy <- 0.28*365 #102.2 days

#up to here

impact_overall <- df_distr_wide_setting %>%
  group_by(ref, init_EIR) %>%
  #group_by(ref) %>%
  filter(between(t, start, start + 365)) %>%
  summarise(tot_cases_baseline = sum(`baseline-Sen`)*1000,
            tot_cases_10d = sum(`10d-stagger-Sen`)*1000,
            tot_cases_20d = sum(`20d-stagger-Sen`)*1000,
            tot_cases_all_in = sum(`all-in-one-Sen`)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_20d = ((tot_cases_baseline - tot_cases_20d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
         scenario = "Once year since start") #measuring overall impact



impact_bohemia <- df_distr_wide_setting %>%
  group_by(ref, init_EIR) %>%
  filter(between(t, start, bohemia_inc_period)) %>%
  summarise(tot_cases_baseline = sum(`baseline-Sen`)*1000,
            tot_cases_10d = sum(`10d-stagger-Sen`)*1000,
            tot_cases_20d = sum(`20d-stagger-Sen`)*1000,
            tot_cases_all_in = sum(`all-in-one-Sen`)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_20d = ((tot_cases_baseline - tot_cases_20d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
         scenario = "bohemia") #measuring overall impact

impact_matamal <- df_distr_wide_setting_prev %>%
  group_by(ref, init_EIR) %>%
  filter(t == matamal_survey) %>%
  summarise(prev_baseline = `baseline-Sen`,
            prev_10d =`10d-stagger-Sen`,
            prev_20d = `20d-stagger-Sen`,
            prev_all_in =`all-in-one-Sen`) %>%
  mutate(impact_all_in = ((prev_baseline - prev_all_in)/prev_baseline)*100,
         impact_20d = ((prev_baseline - prev_20d)/prev_baseline)*100,
         impact_10d = ((prev_baseline - prev_10d)/prev_baseline)*100,
         scenario = "matamal") #measuring overall impact


impact_measurements <- do.call("rbind", list(impact_overall, impact_bohemia, impact_matamal)) %>%
  select(init_EIR, impact_all_in, impact_20d, impact_10d, scenario)

impact_measurements_long <- impact_measurements %>%
  pivot_longer(cols = c(impact_all_in, impact_20d, impact_10d), names_to = "intervention", values_to = "impact")

levels_x <- unique(impact_measurements_long$scenario)

#impact_sens_facet <- ggplot(impact_measurements_long, aes(x = factor(scenario), y = impact, fill = as.factor(intervention)))+
#  geom_bar(stat = "identity", position = position_dodge())+
#  facet_wrap(vars(init_EIR, Q0, ivm_cov), labeller = label_both)+
#  theme_bw()+
#  theme(legend.position = c(0.9, 0.1))+
#  #scale_fill_manual(name = "Time to completel MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
#  #                 labels = c("10 days", "30 days", "1 day"))+
#  xlab("Time period of measurement")+
#  ylab("Difference(%) in endpoint compared to baseline")

#impact for init_EIR = 100 and Q0 = 0.92
distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')


#impact with var in coverage

#impact_plot_data_cov <- impact_measurements_long %>%
#  filter(init_EIR == 100 & Q0 == 0.92 & ivm_cov == covs[2])
#
#impact_plot_error_cov <- impact_measurements_long %>%
#  filter(init_EIR == 100 & Q0 == 0.92 & ivm_cov %in% c(covs[1], covs[3])) %>%
#  mutate(ivm_cov_input = case_when(ivm_cov == covs[1] ~ "impact_low_cov",
#                                   ivm_cov == covs[3] ~ "impact_high_cov")) %>%
#  select(init_EIR, Q0,scenario, intervention, ivm_cov_input, impact) %>%
#  pivot_wider(names_from = ivm_cov_input, values_from = impact)
#
#impact_plot_uncertainty <- left_join(impact_plot_data_cov, impact_plot_error_cov)
#
#
#impact_cov_plot <- ggplot(impact_plot_uncertainty, aes(x = factor(scenario), y = impact, fill = as.factor(intervention)))+
#  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
#  geom_errorbar(aes(ymin = impact_low_cov, ymax = impact_high_cov), position = position_dodge(width = 0.9),
#                width = 0.4)+
#  theme_bw()+
#  scale_fill_manual(name = "Time to completel MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
#                    labels = c("10 days", "20 days", "1 day"))+
#  ylab("Difference (%) in endpoint compared to baseline")+
#  scale_x_discrete(labels = c("bohemia" = "Incidence within \n 6m of start",
#                              "matamal" = "Prevalence 4 weeks \n after last MDA",
#                              "Once year since start" = "Incidence within a \n year since start"),
#                   name = "Time of measurement")+
#  theme(legend.position = c(0.8, 0.9))+
#  ggtitle("Impact with 70% IVM coverage (uncertainty: 50%, 80% cov)")+
#  guides(fill = "none")
#
#
##impact with var in Q0
#
#impact_plot_data_Q0 <- impact_measurements_long %>%
#  filter(init_EIR == 100 & Q0 == 0.92 & ivm_cov == covs[2])
#
#Q0_range <- unique(impact_measurements_long$Q0)
#
#impact_plot_error_Q0 <- impact_measurements_long %>%
#  ungroup() %>%
#  filter(init_EIR == 100 & ivm_cov == covs[2] & Q0 %in% c(Q0_range[1], Q0_range[4])) %>%
#  mutate(Q0_input = case_when(Q0 == Q0_range[1] ~ "impact_low_Q0",
#                              TRUE ~ "impact_high_Q0")) %>%
#  select(init_EIR,scenario, intervention, impact, Q0_input) %>%
#  pivot_wider(names_from = Q0_input, values_from = impact)
#
#impact_plot_uncertainty_Q0 <- left_join(impact_plot_data_Q0, impact_plot_error_Q0)
#
#
#impact_Q0_plot <- ggplot(impact_plot_uncertainty_Q0, aes(x = factor(scenario), y = impact, fill = as.factor(intervention)))+
#  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
#  geom_errorbar(aes(ymin = impact_low_Q0, ymax = impact_high_Q0), position = position_dodge(width = 0.9),
#                width = 0.4)+
#  theme_bw()+
#  scale_fill_manual(name = "Time to completel MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
#                    labels = c("10 days", "20 days", "1 day"))+
#  ylab("Difference (%) in endpoint compared to baseline")+
#  scale_x_discrete(labels = c("bohemia" = "Incidence within \n 6m of start",
#                              "matamal" = "Prevalence 4 weeks \n after last MDA",
#                              "Once year since start" = "Incidence within a \n year since start"),
#                   name = "Time of measurement")+
#  theme(legend.position = c(0.8, 0.9))+
#  ggtitle("Impact with 0.92 Q0 (uncertainty: 0.21, 0.94)")+
#  guides(fill = "none")+
#  ylim(0,60)




impact_main_plot <- impact_measurements_long %>%
  filter(init_EIR == 100) %>%
  ggplot()+
  aes(x = factor(scenario), y = impact, fill = as.factor(intervention))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_bw()+
  theme(legend.position = c(0.9, 0.1))+
  scale_fill_manual(name = "Time to completel MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
                    labels = c("10 days", "30 days", "1 day"))+
  xlab("Time of measurement")+
  ylab("Difference(%) in endpoint compared to baseline")+
  #ylim(0, 60)+
  guides(fill = "none")+
  #facet_wrap(vars(ivm_cov))+
  scale_x_discrete(labels = c("bohemia" = "Incidence U5s (start to 6m later)",
                              "matamal" = "Prevalence U5s (1m after last MDA)",
                              "Once year since start" = "Incidence U5s (start to 1y later)"))


distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
dynamics_Sen <- cowplot::plot_grid(mv_plot, eir_plot, prev_plot, inc_plot,
                                         labels = c("A", "B", "C", "D"),
                                         align = "v")

#impact_perennial <- cowplot::plot_grid(impact_cov_plot, impact_Q0_plot, labels = c("E", "F"),
#                                      nrow = 2, align = "v")


plot_seasonal<- cowplot::plot_grid(dynamics_Sen, impact_main_plot,
                                   labels = c("", "E"))

#plot_perennial_HS <- cowplot::plot_grid(dynamics_perennial, impact_cov_plot,
#                                        labels = c("","E"))

#ggsave(plot_perennial, file = "analysis/target-profiles-distrib-strat/chapter/plots/impact_perennial.svg")
#ggsave(plot_perennial, file = "analysis/target-profiles-distrib-strat/chapter/plots/impact_perennial.pdf")

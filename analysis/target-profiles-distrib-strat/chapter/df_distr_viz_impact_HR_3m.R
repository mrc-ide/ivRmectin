#viz and epi impact of staggered vs all in one distribution

#sensitivity analysis: endemicity, (init_EIR = 2 or 100), Q0 (four levels, for each species) and seasonality (no seasonality or Fattick, Senegal)
require(tidyverse)

start <- (365*5)+200

df_distr_all <- readRDS("analysis/target-profiles-distrib-strat/chapter/output/df_distr_HR_3m.rds")

df_distr_all %>%
  filter(t == 1) %>% #EIR at t = 1 either 1.39 or 69.5
  group_by(ref) %>%
  summarise(eir = EIR_tot) %>%
  distinct() #so ref 1-4 are low EIR and ref 5-8 are high EIR

df_distr_all <- df_distr_all %>%
  mutate(init_EIR = case_when(ref <= 4 ~ 2,
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

df_distr <- df_distr_all %>%
  filter(init_EIR == 100 & Q0 == 0.92)


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
  ylim(0, 50) +
  ylab("Mosquito density")+
  theme(legend.position = c(0.6, 0.2)) +

  scale_color_manual(name = "Scenario", values = distr_pals,
                     labels = c("10 days to complete 1 round of MDA", "30 days to complete 1 round of MDA",
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

mosq_killed <- ggplot(df_distr, aes(x = (t - start)/365, y = mvx_dead, col = as.factor(model_type))) +
  geom_line() +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +

  scale_color_manual(name = "Scenario", values = distr_pals) +
  xlim(-0.25, 1)+
  guides(col = "none")+
  ylab("Number of mosquitoes killed by intervention")+
  xlab("Years since intervention started")+
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 1, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = plot_matamal/365, col = "grey", linetype = "dashed", size = 1.05)  +
  geom_vline(xintercept = plot_bohemia/365, col = "orange", linetype = "dashed", size = 1.05)

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
  ylim(0, 100)+
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
  ylim(0, 75)+
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
  ylim(0, 7.5)+
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

dynamics <- cowplot::plot_grid(mv_plot, mosq_killed, eir_plot, inc_plot, prev_plot, align = "v")

#then look at epi impact
#df_int <- do.call("rbind", list(df_mod1, df_mod2, df_mod3))

#by day 300, all back to eqm for mv

#look at epi impact for all

#need to look at this over different time periods. Note, that this is over 5y so overall impact is very diluted.
#slower strategy takes longer to return to eqm hence appears more impactful

#make df_distr_all according to each scenario

df_distr_wide_setting <- df_distr_all %>%
  select(t, Q0, init_EIR, model_type, clin_inc0to5) %>%
  group_by(Q0, init_EIR) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)
  })


names(df_distr_wide_setting)

df_distr_wide_setting_prev <- df_distr_all %>%
  select(t, Q0, init_EIR, model_type, slide_prev0to5) %>%
  group_by(Q0, init_EIR) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to5)
  })

#measure impact from start to 1y later
#timy <- 0.28*365 #102.2 days

impact_overall <- df_distr_wide_setting %>%
  group_by(init_EIR, Q0) %>%
  filter(between(t, start, start + 365)) %>%
  summarise(tot_cases_baseline = sum(baseline)*1000,
            tot_cases_10d = sum(`10d-stagger`)*1000,
            tot_cases_30d = sum(`30d-stagger`)*1000,
            tot_cases_all_in = sum(`all-in-one`)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
         scenario = "Once year since start") #measuring overall impact


impact_bohemia <- df_distr_wide_setting %>%
  group_by(init_EIR, Q0) %>%
  filter(between(t, start, bohemia_inc_period)) %>%
  summarise(tot_cases_baseline = sum(baseline)*1000,
            tot_cases_10d = sum(`10d-stagger`)*1000,
            tot_cases_30d = sum(`30d-stagger`)*1000,
            tot_cases_all_in = sum(`all-in-one`)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
         scenario = "bohemia") #measuring overall impact

impact_matamal <- df_distr_wide_setting_prev %>%
  group_by(init_EIR, Q0) %>%
  filter(t == matamal_survey) %>%
  summarise(prev_baseline = baseline,
            prev_10d =`10d-stagger`,
            prev_30d = `30d-stagger`,
            prev_all_in =`all-in-one`) %>%
  mutate(impact_all_in = ((prev_baseline - prev_all_in)/prev_baseline)*100,
         impact_30d = ((prev_baseline - prev_30d)/prev_baseline)*100,
         impact_10d = ((prev_baseline - prev_10d)/prev_baseline)*100,
         scenario = "matamal") #measuring overall impact


impact_measurements <- do.call("rbind", list(impact_overall, impact_bohemia, impact_matamal)) %>%
  select(init_EIR, Q0, impact_all_in, impact_30d, impact_10d, scenario)

impact_measurements_long <- impact_measurements %>%
  pivot_longer(cols = c(impact_all_in, impact_30d, impact_10d), names_to = "intervention", values_to = "impact")

levels_x <- unique(impact_measurements_long$scenario)

impact_sens_facet <- ggplot(impact_measurements_long, aes(x = factor(scenario), y = impact, fill = as.factor(intervention)))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(init_EIR, Q0), labeller = label_both)+
  theme_bw()+
  theme(legend.position = c(0.9, 0.1))+
  #scale_fill_manual(name = "Time to completel MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
   #                 labels = c("10 days", "30 days", "1 day"))+
  xlab("Time period of measurement")+
  ylab("Difference(%) in endpoint compared to baseline")

#impact for init_EIR = 100 and Q0 = 0.92
distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')

impact_main_plot <- impact_measurements_long %>%
  filter(init_EIR == 100 & Q0 == 0.92) %>%
  ggplot()+
  aes(x = factor(scenario), y = impact, fill = as.factor(intervention))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_bw()+
  theme(legend.position = c(0.9, 0.1))+
  scale_fill_manual(name = "Time to completel MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
                    labels = c("10 days", "30 days", "1 day"))+
  xlab("Time of measurement")+
  ylab("Difference(%) in endpoint compared to baseline")+
  ylim(0, 60)+
  guides(fill = "none")+
  scale_x_discrete(labels = c("bohemia" = "BOHEMIA protocol",
                              "matamal" = "MATAMAL protocol",
                              "Once year since start" = "One year since start"))

distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
dynamics_perennial <- cowplot::plot_grid(mv_plot, mosq_killed, eir_plot, prev_plot, inc_plot,
                                         impact_main_plot, labels = c("A", "B", "C", "D", "E", "F"),
                                         align = "v")



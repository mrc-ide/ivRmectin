#viz and epi impact of staggered vs all in one distribution

#sensitivity analysis: endemicity, (init_EIR = 2 or 100), Q0 (four levels, for each species) and seasonality (no seasonality or Fattick, Senegal)

start <- (365*5)+200

df_distr_all <- readRDS("analysis/target-profiles-distrib-strat/chapter/output/df_distr_HR_Sen.rds")

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
  ylim(0, 200) +
  ylab("Mosquito density")+
  theme(legend.position = c(0.7, 0.5)) +

  scale_color_manual(name = "Scenario", values = distr_pals) +
  xlim(-0.25, 1) +
  xlab("Years since intervention started") +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 0.28, col = "black", linetype = "dashed", size = 1.05)  #+

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
  geom_line(size = 1.1) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +

  scale_color_manual(name = "Scenario", values = distr_pals) +
  xlim(-0.25, 1)+
  guides(col = "none")+
  ylab("Number of mosquitoes killed by intervention")+
  xlab("Years since intervention started")+
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 0.28, col = "black", linetype = "dashed", size = 1.05) #+

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
  #ylim(0, 100)+
  scale_color_manual(name = "Scenario", values = distr_pals)+
  guides(col = "none")+
  xlim(-0.25, 1) +
  ylab("Average annual EIR")+
  xlab("Years since intervention started") +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 0.28, col = "black", linetype = "dashed", size = 1.05) #+

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
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)+
  guides(col = "none")+
  xlim(-0.25, 1) +
  xlab("Years since intervention started")+
  ylab("Slide prevalence (%) in children \n under 5-years-old")+
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 0.28, col = "black", linetype = "dashed", size = 1.05) #+

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
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)+
  guides(col = "none")+
  ylab("Clinical incidence in children \n under 5-years-old, per 1000 persons")+
  #ylim(0, 7.5)+
  xlim(-0.25, 1) +
  xlab("Years since intervention started")+
  geom_vline(xintercept = 0, col = "black", linetype = "dashed", size = 1.05) +
  geom_vline(xintercept = 0.28, col = "black", linetype = "dashed", size = 1.05) #+

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

dynamics <- cowplot::plot_grid(mv_plot, mosq_killed, eir_plot, inc_plot, prev_plot)

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

#measure impact between t start and returns to eqm 0.28y since int started
timy <- 0.28*365 #102.2 days

impact_overall <- df_distr_wide_setting %>%
  group_by(init_EIR, Q0) %>%
  filter(between(t, start, start + timy)) %>%
  summarise(tot_cases_baseline = sum(`baseline-Sen`)*1000,
            tot_cases_10d = sum(`10d-stagger-Sen`)*1000,
            tot_cases_30d = sum(`30d-stagger-Sen`)*1000,
            tot_cases_all_in = sum(`all-in-one-Sen`)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
         scenario = "overall") #measuring overall impact


#impact in different time periods. 10d and 30d - all distr are done and all mosquitoes killed by both but non-linearity and lags
#means that the clinical impact of this is not realised yet so 10d intervention looks better
impact_30_stag <- df_distr_wide_setting %>%
  group_by(init_EIR, Q0) %>%
  filter(between(t, start, start+53)) %>% #in population for 53 days
  summarise(tot_cases_baseline = sum(`baseline-Sen`)*1000,
            tot_cases_10d = sum(`10d-stagger-Sen`)*1000,
            tot_cases_30d = sum(`30d-stagger-Sen`)*1000,
            tot_cases_all_in = sum(`all-in-one-Sen`)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
         scenario = "30-day-measurement")

#not completed the MDA for the 30-day staggered delivery yet, so 10-day appears better
impact_10_stag <- df_distr_wide_setting %>%
  group_by(init_EIR, Q0) %>%
  filter(between(t, start, start+33)) %>% #in pop for 33 days
  summarise(tot_cases_baseline = sum(`baseline-Sen`)*1000,
            tot_cases_10d = sum(`10d-stagger-Sen`)*1000,
            tot_cases_30d = sum(`30d-stagger-Sen`)*1000,
            tot_cases_all_in = sum(`all-in-one-Sen`)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
         scenario = "10-day-measurement")

impact_measurements <- do.call("rbind", list(impact_overall, impact_30_stag, impact_10_stag)) %>%
  select(init_EIR, Q0, impact_all_in, impact_30d, impact_10d, scenario)

impact_measurements_long <- impact_measurements %>%
  pivot_longer(cols = c(impact_all_in, impact_30d, impact_10d), names_to = "intervention", values_to = "cases_averted")

levels_x <- unique(impact_measurements_long$scenario)

impact_sens_facet <- ggplot(impact_measurements_long, aes(x = factor(scenario, levels = c(levels_x[1], levels_x[3], levels_x[2]),
                                                                     labels = c("Overall", "33d since start", "53d since start")), y = cases_averted, fill = as.factor(intervention)))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(init_EIR, Q0), labeller = label_both)+
  theme_bw()+
  theme(legend.position = c(0.9, 0.1))+
  scale_fill_manual(name = "Time to complete MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
                    labels = c("10 days", "30 days", "1 day"))+
  xlab("Time period of measurement")+
  ylab("Clinical cases averted (%) in children under 5-years-old")

#impact for init_EIR = 100 and Q0 = 0.92
impact_main_plot <- impact_measurements_long %>%
  filter(init_EIR == 100 & Q0 == 0.92) %>%
  ggplot()+
  aes(x = factor(scenario, levels = c(levels_x[1], levels_x[3], levels_x[2]),
                 labels = c("Overall", "33d since start", "53d since int")), y = cases_averted, fill = as.factor(intervention))+
  geom_bar(stat = "identity", position = position_dodge())+
  theme_bw()+
  theme(legend.position = c(0.9, 0.1))+
  scale_fill_manual(name = "Time to completel MDA", values = c(distr_pals[1], distr_pals[2], distr_pals[3]),
                    labels = c("10 days", "30 days", "1 day"))+
  xlab("Time period of measurement")+
  ylab("Clinical cases averted (%) in \n children under 5-years-old")+
  guides(fill = "none")

distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
dynamics_seasonal <- cowplot::plot_grid(mv_plot, mosq_killed, eir_plot, prev_plot, inc_plot,
                                         impact_main_plot, labels = c("A", "B", "C", "D", "E", "F"),
                                         align = "v")


#then for prevalence, but look at set time points
df_distr_wide_prev_setting <- df_distr_all %>%
  select(t, Q0, init_EIR, model_type, slide_prev0to5) %>%
  group_by(Q0, init_EIR) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to5)
  })

#take prevalence estimate at 3 months after MDA started

impact_end_prev <- df_distr_wide_prev_setting %>%
  group_by(init_EIR, Q0) %>%
  filter(t == start + (3*30)) %>% #take estimate 3 months after MDA started
  summarise(prev_baseline = baseline,
            prev_10d = `stag-10d`,
            prev_30d = `stag-30d`,
            prev_all_in = `all-in`) %>%
  mutate(impact_all_in = ((prev_baseline - prev_all_in)/prev_baseline)*100,
         impact_30d = ((prev_baseline - prev_30d)/prev_baseline)*100,
         impact_10d = ((prev_baseline - prev_10d)/prev_baseline)*100,
         scenario = "point-prev-3months-post-MDA") #measuring overall impact

#plot number of mosquitoes killed over time - the higher coverage of the non-staggered method causes a greater total number of mosq to be killed
#diluting the coverage means fewer mosquitoes are killed by the staggered approach, but the difference in amount of time to kill that target number
#does not lead to detectable differences in impact.


#viz and epi impact of staggered vs all in one distribution

#sensitivity analysis: endemicity, (init_EIR = 2 or 100), Q0 (four levels, for each species) and seasonality (no seasonality or Fattick, Senegal)

start <- 200

df_distr_all <- readRDS("analysis/target-profiles-distrib-strat/chapter/output/df_distr.rds")

df_distr_all %>%
  filter(t == 1) %>% #EIR at t = 1 either 1.39 or 69.5
  group_by(ref) %>%
  summarise(eir = EIR_tot) %>%
  distinct() #so ref 1-4 are low EIR and ref 5-8 are high EIR

df_distr_all <- df_distr_all %>%
  mutate(init_EIR = case_when(ref <= 4 ~ 2,
                              TRUE ~ 100))

unique(df_distr_all$init_EIR)


distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')

#filter to just high EIR and gamb-like vector (Q0 = 0.92) for plots

df_distr <- df_distr_all %>%
  filter(init_EIR == 100 & Q0 == 0.92)

#check pop has been dividided up correctly, can see stagger

ggplot(df_distr, aes(x = t, y = mvtot_1, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  theme(legend.position = c(0.7, 0.2))+
  ylim(0, 15)+
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)

ggplot(df_distr, aes(x = t, y = mvtot_2, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  ylim(0, 15)+
  theme(legend.position = c(0.7, 0.2))+
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)

ggplot(df_distr, aes(x = t, y = mvtot_3, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  ylim(0, 15)+
  theme(legend.position = c(0.7, 0.2))+
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)


mv_plot <- ggplot(df_distr, aes(x = t, y = mv, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  ylim(0, 50)+
  theme(legend.position = c(0.7, 0.2))+
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)


eir_plot <- ggplot(df_distr, aes(x = t, y = EIR_tot, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  ylim(0, 100)+
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)+
  guides(col = "none")

prev_plot <- ggplot(df_distr, aes(x = t, y = slide_prev0to5, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  ylim(0, 1)+
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)+
  guides(col = "none")

inc_plot <- ggplot(df_distr, aes(x = t, y = clin_inc0to5*1000, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw()+
  scale_color_manual(name = "Time to complete MDA", values = distr_pals)+
  guides(col = "none")+
  ylab("clinical incidence u5 per 1000 persons")+
  ylim(0, 10)


cowplot::plot_grid(mv_plot, eir_plot, prev_plot, inc_plot, labels = c("A", "B", "C", "D"))

#then look at epi impact
#df_int <- do.call("rbind", list(df_mod1, df_mod2, df_mod3))

#by day 300, all back to eqm for mv

#look at epi impact for all

df_distr_wide <- df_distr %>%
  filter(init_EIR == 100 & Q0 == 0.92) %>%
  select(t, model_type, clin_inc0to5) %>%
  pivot_wider(names_from = "model_type", values_from = clin_inc0to5)

#need to look at this over different time periods. Note, that this is over 5y so overall impact is very diluted.
#slower strategy takes longer to return to eqm hence appears more impactful
impact <- df_distr_wide %>%
  summarise(tot_cases_baseline = sum(baseline)*1000,
            tot_cases_10d = sum(`stag-10d`, na.rm = TRUE)*1000,
            tot_cases_30d = sum(`stag-30d`, na.rm = TRUE)*1000,
            tot_cases_all_in = sum(`all-in`, na.rm = TRUE)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100)


#more automated version
impact_setting <- df_distr_all %>%
  select(t, Q0, init_EIR, model_type, clin_inc0to5) %>%
  group_by(Q0, init_EIR) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)

    tibble(
      tot_cases_baseline = sum(df_wide$baseline, na.rm = TRUE) * 1000,
      tot_cases_10d = sum(df_wide$`stag-10d`, na.rm = TRUE) * 1000,
      tot_cases_30d = sum(df_wide$`stag-30d`, na.rm = TRUE) * 1000,
      tot_cases_all_in = sum(df_wide$`all-in`, na.rm = TRUE) * 1000
    ) %>%
      mutate(
        impact_all_in = ((tot_cases_baseline - tot_cases_all_in) / tot_cases_baseline) * 100,
        impact_30d = ((tot_cases_baseline - tot_cases_30d) / tot_cases_baseline) * 100,
        impact_10d = ((tot_cases_baseline - tot_cases_10d) / tot_cases_baseline) * 100
      )
  }) %>%
  ungroup()

#impact in different time periods. 10d spacing - all distr are done, but not for the 30d spacing, the impact of last one not seen yet
#so 10d appears to be better
impact_30d <- df_distr_wide %>%
  filter(between(t, start, start+30)) %>%
  summarise(tot_cases_baseline = sum(baseline, na.rm = TRUE)*1000,
            tot_cases_10d = sum(`stag-10d`, na.rm = TRUE)*1000,
            tot_cases_30d = sum(`stag-30d`, na.rm = TRUE)*1000,
            tot_cases_all_in = sum(`all-in`, na.rm = TRUE)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100)

#not seeing that the green line is returning to eqm quicker at this time point, and drop in green has been greater so green appears bigger
impact_55d <- df_distr_wide %>%
  filter(between(t, start, start+55)) %>%
  summarise(tot_cases_baseline = sum(baseline, na.rm = TRUE)*1000,
            tot_cases_10d = sum(`stag-10d`, na.rm = TRUE)*1000,
            tot_cases_30d = sum(`stag-30d`, na.rm = TRUE)*1000,
            tot_cases_all_in = sum(`all-in`, na.rm = TRUE)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100)

#still see impact of reductions for the staggered distribution
impact_200d <- df_distr_wide %>%
  filter(between(t, start, start+100)) %>%
  summarise(tot_cases_baseline = sum(baseline, na.rm = TRUE)*1000,
            tot_cases_10d = sum(`stag-10d`, na.rm = TRUE)*1000,
            tot_cases_30d = sum(`stag-30d`, na.rm = TRUE)*1000,
            tot_cases_all_in = sum(`all-in`, na.rm = TRUE)*1000) %>%
  mutate(impact_all_in = ((tot_cases_baseline - tot_cases_all_in)/tot_cases_baseline)*100,
         impact_30d = ((tot_cases_baseline - tot_cases_30d)/tot_cases_baseline)*100,
         impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100)

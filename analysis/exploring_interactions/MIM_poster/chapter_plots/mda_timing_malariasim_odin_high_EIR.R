require(tidyverse)

#working out the correct endec_mu and wane values for malariasim_exp_decay.R

#read in odin models
odin_mda_early <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_early.rds")
odin_mda_medium <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_medium.rds")
odin_mda_late <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_late.rds")
odin_baseline <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_baseline.rds")
odin_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/chapter_plots/df_ITN.rds")

#read in malariasim exp decay models
msim_mda_early <- readRDS("W:/endectocides-cluster/data/df_var1_mda_early_b_malariasim_odin.rds")
msim_mda_medium <- readRDS("W:/endectocides-cluster/data/df_var1_mda_medium_b_malariasim_odin.rds")
msim_mda_late <- readRDS("W:/endectocides-cluster/data/df_var1_mda_late_b_malariasim_odin.rds")
msim_baseline <- readRDS("W:/endectocides-cluster/data/df_baselineb_malariasim_odin.rds")
msim_ITN <- readRDS("W:/endectocides-cluster/data/df_df_var1_ITNb_malariasim_odin.rds")

#use least-squares fitting to calculate wane and endec_mu between models
df_var1 <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-mda-times.rds")
high_EIR <- unique(df_var1$init_EIR)[2]
#first filter in period of ivermectin MDA
time_period <- 365*15
mda_int <- 30

#introduce nets 1y into simulation

itn_on <- 365*5

net_seq <- seq(itn_on, time_period, by = 3*365)


#IVM MDA starts 6m after last ITN campaign
IVM_begin1 <- net_seq[2]+180
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

#IVM MDA starts 1y after last ITN campaign
IVM_begin2 <- net_seq[2] + 365
IVM_start2 <- c(IVM_begin2, IVM_begin2+mda_int, IVM_begin2+mda_int+mda_int)

#IVM starts 2.5 after last ITN campaign
IVM_begin3 <- net_seq[2] + (2*365)
IVM_start3 <- c(IVM_begin3, IVM_begin3+mda_int, IVM_begin3+mda_int+mda_int)

#parameter set for endec_mu and wane
endec_mu_vec <- seq(0.05,0.15,0.01)
wane_vec <- seq(0,0.1,0.01)
combos <- expand.grid(endec_mu_vec, wane_vec)


#early####
early_t_start <- IVM_start1[1]
early_t_end <-  IVM_start1[3]+23

# ref == 2 is high EIR

time_odin_mda_early <- odin_mda_early %>%
  filter(ref == 2 & between(t, early_t_start, early_t_end)) %>%
  select(t,Ivtot)

check_ref <- msim_mda_early %>%
  group_by(ref) %>%
  filter(t == 1) %>%
  summarise(prev = slide_prev0to5)

#ref 1 to 121 is low EIR
#ref 122 to 242 is high EIR

time_msim_mda_early <- msim_mda_early %>%
  filter(ref >= 122 & ref <= 242 & between(t, early_t_start, early_t_end)) %>%
  select(t,Ivtot, ref)
unique(time_msim_mda_early$t)

list_msim_early <- split(time_msim_mda_early, f = time_msim_mda_early$ref)

error_early <- numeric()
for (i in 1:nrow(combos)){
  error_early <- c(error_early, sum((time_odin_mda_early$Ivtot - list_msim_early[[i]]$Ivtot)^2))
}
range(error_early)
index_early <- which.min(error_early)
best_fit_early <- combos[index_early,]

#best fit is endec_mu of 0.14 and wane of 0.01 Need to double check that this cross-refs correctly to the endec_mu and wane in that model

#medium distr###

medium_t_start <- IVM_start2[1]
medium_t_end <-  IVM_start2[3]+23

# ref == 2 is high EIR

time_odin_mda_medium <- odin_mda_medium %>%
  filter(ref == 2 & between(t, medium_t_start, medium_t_end)) %>%
  select(t,Ivtot)

check_ref <- msim_mda_medium %>%
  group_by(ref) %>%
  filter(t == 1) %>%
  summarise(prev = slide_prev0to5)

#ref 1 to 121 is low EIR
#ref 122 to 242 is high EIR

time_msim_mda_medium <- msim_mda_medium %>%
  filter(ref >= 122 & ref <= 242 & between(t, medium_t_start, medium_t_end)) %>%
  select(t,Ivtot, ref)


list_msim_medium <- split(time_msim_mda_medium, f = time_msim_mda_medium$ref)

error_medium <- numeric()
for (i in 1:nrow(combos)){
  error_medium <- c(error_medium, sum((time_odin_mda_medium$Ivtot - list_msim_medium[[i]]$Ivtot)^2))
}
range(error_medium)
index_medium <- which.min(error_medium)
best_fit_medium <- combos[index_medium,]
#also endec_mu = 0.14 and wane = 0.01

#late####
late_t_start <- IVM_start3[1]
late_t_end <-  IVM_start3[3]+23

# ref == 2 is high EIR

time_odin_mda_late <- odin_mda_late %>%
  filter(ref == 2 & between(t, late_t_start, late_t_end)) %>%
  select(t,Ivtot)

check_ref <- msim_mda_late %>%
  group_by(ref) %>%
  filter(t == 1) %>%
  summarise(prev = slide_prev0to5)

#ref 1 to 121 is low EIR
#ref 122 to 242 is high EIR

time_msim_mda_late <- msim_mda_late %>%
  filter(ref >= 122 & ref <= 242 & between(t, late_t_start, late_t_end)) %>%
  select(t,Ivtot, ref)

list_msim_late <- split(time_msim_mda_late, f = time_msim_mda_late$ref)

error_late <- numeric()
for (i in 1:nrow(combos)){
  error_late <- c(error_late, sum((time_odin_mda_late$Ivtot - list_msim_late[[i]]$Ivtot)^2))
}
range(error_late)
index_late <- which.min(error_late)
best_fit_late <- combos[index_late,]

#endec mu is 0.14 for all and wane is 0.01
#could restrict search so from endec_mu 0.05 to 0.15 and wane 0 to 0.01.

#filter for each scenario
endec_mu_vec

msim_early_bestfit <- msim_mda_early %>%
  filter(ref >= 122 & endec_mu == endec_mu_vec[10] & wane == wane_vec[2])

msim_medium_bestfit <- msim_mda_medium %>%
  filter(ref >= 122&endec_mu == endec_mu_vec[10] & wane == wane_vec[2])

msim_late_bestfit <- msim_mda_late %>%
  filter(ref >= 122&endec_mu == endec_mu_vec[10] & wane == wane_vec[2])

high_eir_odin_mda_early <- odin_mda_early %>%
  filter(ref == 2)

high_eir_odin_mda_medium <- odin_mda_medium %>%
  filter(ref == 2)

high_eir_odin_mda_late <- odin_mda_late %>%
  filter(ref == 2)

vars <- c("Ivtot", "EIRout", "slide_prev0to5", "clin_inc0to5")

mod_pals <-  c('#1b9e77', '#d95f02')

#early
plots_early <- lapply(vars, function(v) {
  if (v == "EIRout") {
    p <- ggplot(high_eir_odin_mda_early, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000)) +
      geom_line(col =mod_pals[1]) +
      geom_line(data = msim_early_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000), col = mod_pals[2]) +
      labs(y = "Daily EIR (per 1000 persons)", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else if (v == "clin_inc0to5") {
    p <- ggplot(high_eir_odin_mda_early, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000)) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_early_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000), col = mod_pals[2]) +
      labs(y = "Clinical incidence \n in under 5-year-olds (per 1000 persons)", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else if (v == "slide_prev0to5") {
    p <- ggplot(high_eir_odin_mda_early, aes(x = (t-net_seq[1])/365, y = .data[[v]]*100)) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_early_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]*100), col = mod_pals[2]) +
      labs(y = "Slide prevalence (%) in \n under 5-year-olds", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5), ylim = c(0,80))

  } else if (v == "Ivtot") {
    p <- ggplot(high_eir_odin_mda_early, aes(x = (t-net_seq[1])/365, y = .data[[v]])) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_early_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]), col = mod_pals[2]) +
      labs(y = "Infectious vectors", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else {
    # fallback (shouldn’t really happen here)
    p <- ggplot(high_eir_odin_mda_early, aes(x = (t-net_seq[1])/365, y = .data[[v]])) +
      geom_line() +
      geom_line(data = msim_early_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]), col = "red") +
      labs(title = v) +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))
  }
  p
})

early_fits <- cowplot::plot_grid(plotlist = plots_early,
                                labels = c("A", "B", "C", "D"))

ggsave(early_fits,file = "analysis/exploring_interactions/MIM_poster/chapter_plots/early_fits_odin_msimodin_high_eir.pdf")

#medium

plots_medium <- lapply(vars, function(v) {
  if (v == "EIRout") {
    p <- ggplot(high_eir_odin_mda_medium, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000)) +
      geom_line(col =mod_pals[1]) +
      geom_line(data = msim_medium_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000), col = mod_pals[2]) +
      labs(y = "Daily EIR (per 1000 persons)", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else if (v == "clin_inc0to5") {
    p <- ggplot(high_eir_odin_mda_medium, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000)) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_medium_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000), col = mod_pals[2]) +
      labs(y = "Clinical incidence \n in under 5-year-olds (per 1000 persons)", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else if (v == "slide_prev0to5") {
    p <- ggplot(high_eir_odin_mda_medium, aes(x = (t-net_seq[1])/365, y = .data[[v]]*100)) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_medium_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]*100), col = mod_pals[2]) +
      labs(y = "Slide prevalence (%) in \n under 5-year-olds", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5), ylim = c(0,80))

  } else if (v == "Ivtot") {
    p <- ggplot(high_eir_odin_mda_medium, aes(x = (t-net_seq[1])/365, y = .data[[v]])) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_medium_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]), col = mod_pals[2]) +
      labs(y = "Infectious vectors", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else {
    # fallback (shouldn’t really happen here)
    p <- ggplot(high_eir_odin_mda_medium, aes(x = (t-net_seq[1])/365, y = .data[[v]])) +
      geom_line() +
      geom_line(data = msim_medium_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]), col = "red") +
      labs(title = v) +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))
  }
  p
})

medium_fits <- cowplot::plot_grid(plotlist = plots_medium,
                                 labels = c("A", "B", "C", "D"))

ggsave(medium_fits,file = "analysis/exploring_interactions/MIM_poster/chapter_plots/medium_fits_odin_msimodin_high_eir.pdf")


#late


plots_late <- lapply(vars, function(v) {
  if (v == "EIRout") {
    p <- ggplot(high_eir_odin_mda_late, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000)) +
      geom_line(col =mod_pals[1]) +
      geom_line(data = msim_late_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000), col = mod_pals[2]) +
      labs(y = "Daily EIR (per 1000 persons)", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else if (v == "clin_inc0to5") {
    p <- ggplot(high_eir_odin_mda_late, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000)) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_late_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]] * 1000), col = mod_pals[2]) +
      labs(y = "Clinical incidence \n in under 5-year-olds (per 1000 persons)", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else if (v == "slide_prev0to5") {
    p <- ggplot(high_eir_odin_mda_late, aes(x = (t-net_seq[1])/365, y = .data[[v]]*100)) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_late_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]*100), col = mod_pals[2]) +
      labs(y = "Slide prevalence (%) in \n under 5-year-olds", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5), ylim = c(0,80))

  } else if (v == "Ivtot") {
    p <- ggplot(high_eir_odin_mda_late, aes(x = (t-net_seq[1])/365, y = .data[[v]])) +
      geom_line(col = mod_pals[1]) +
      geom_line(data = msim_late_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]), col = mod_pals[2]) +
      labs(y = "Infectious vectors", x = "Time since first ITN campaign (years)") +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))

  } else {
    # fallback (shouldn’t really happen here)
    p <- ggplot(high_eir_odin_mda_late, aes(x = (t-net_seq[1])/365, y = .data[[v]])) +
      geom_line() +
      geom_line(data = msim_late_bestfit, aes(x = (t-net_seq[1])/365, y = .data[[v]]), col = "red") +
      labs(title = v) +
      theme_bw() +
      coord_cartesian(xlim = c(-0.5, 7.5))
  }
  p
})

late_fits <- cowplot::plot_grid(plotlist = plots_late,
                                labels = c("A", "B", "C", "D"))

ggsave(late_fits,file = "analysis/exploring_interactions/MIM_poster/chapter_plots/late_fits_odin_msimodin_high_eir.pdf")

ggplot(high_eir_odin_mda_late, aes(x = t, y = clin_inc0to5 * 1000)) +
  geom_line() +
  geom_line(data = msim_late_bestfit, aes(x = t, y = clin_inc0to5 * 1000), col = "red") +
  theme_bw()+
  coord_cartesian(xlim = c(3000, 4000))

#compute difference in efficacy (Ivtot, prevalence and clinical incidence for each of these)
#for Ivtot, in the killing period
#prevalence: 1 month after last MDA
#clinical incidence: from start of MDa to 6m later

#filter time appropriately.
high_eir_msim_baseline <- msim_baseline %>%
  filter(ref ==2)

ggplot(high_eir_msim_baseline, aes(x = t, y = Ivtot))+
  geom_line()



baseline_Ivtot <-high_eir_msim_baseline %>%
  filter(between(t, late_t_start, late_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_baseline = Ivtot)

late_mda_Ivtot <- msim_late_bestfit %>%
  filter(between(t, late_t_start, late_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_int = Ivtot)

late_mda_comp <- left_join(baseline_Ivtot, late_mda_Ivtot)

late_mda_eff <- late_mda_comp %>%
  summarise(diff_inf_vec = ((sum(Ivtot_baseline) - sum(Ivtot_int))/ sum(Ivtot_baseline))*100)

#compared to odin model
high_eir_odin <- odin_baseline %>%
  filter(ref == 2 & between(t, late_t_start, late_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_baseline = Ivtot)

high_eir_odin_late <- odin_mda_late %>%
  filter(ref ==2 & between(t, late_t_start, late_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_int = Ivtot)

late_mda_comp_odin <- left_join(high_eir_odin, high_eir_odin_late)

late_mda_odin_eff <- late_mda_comp_odin %>%
  summarise(diff_inf_vec = ((sum(Ivtot_baseline) - sum(Ivtot_int))/ sum(Ivtot_baseline))*100)

#eff diff in prevalence
late_mda_prev <- msim_late_bestfit %>%
  filter(t == late_t_end+30&  endec_mu == 0.14 & wane == 0.01 ) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_int = slide_prev0to5)

baseline_prev <-high_eir_msim_baseline %>%
  filter(t == late_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_baseline = slide_prev0to5)

late_mda_comp <- left_join(baseline_prev, late_mda_prev)
late_mda_eff_prev <- late_mda_comp %>%
  summarise(diff_prev = ((sum(prev_baseline) - sum(prev_int))/ sum(prev_baseline))*100)

#compared to odin model
high_eir_odin_prev <- odin_baseline %>%
  filter(ref == 2 & t == late_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_baseline = slide_prev0to5)

high_eir_odin_late_prev <- odin_mda_late %>%
  filter(ref ==2 & t == late_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_int = slide_prev0to5)

late_mda_comp_odin <- left_join(high_eir_odin_prev, high_eir_odin_late_prev)

late_mda_odin_eff_prev <- late_mda_comp_odin %>%
  summarise(diff_prev = ((sum(prev_baseline) - sum(prev_int))/ sum(prev_baseline))*100)

#eff diff inc u5
late_mda_inc <- msim_late_bestfit %>%
  filter(between(t, late_t_start, late_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_int = clin_inc0to5)

baseline_inc <- high_eir_msim_baseline %>%
  filter(between(t, late_t_start, late_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_baseline = clin_inc0to5)

late_mda_comp <- left_join(baseline_inc, late_mda_inc)
late_mda_eff_inc <- late_mda_comp %>%
  summarise(diff_inc = ((sum(inc_baseline) - sum(inc_int))/ sum(inc_baseline))*100)

#compared to odin model
high_eir_odin_inc <- odin_baseline %>%
  filter(ref == 2 & between(t, late_t_start, late_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_baseline = clin_inc0to5)

high_eir_odin_late_inc <- odin_mda_late %>%
  filter(ref ==2 & between(t, late_t_start, late_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_int = clin_inc0to5)

late_mda_comp_odin <- left_join(high_eir_odin_inc, high_eir_odin_late_inc)

late_mda_odin_eff_inc <- late_mda_comp_odin %>%
  summarise(diff_inc = ((sum(inc_baseline) - sum(inc_int))/ sum(inc_baseline))*100)

#repeat for medium
medium_mda_Ivtot <- msim_medium_bestfit %>%
  filter(between(t, medium_t_start, medium_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_int = Ivtot)

medium_baseline_Ivtot <-high_eir_msim_baseline %>%
  filter(between(t, medium_t_start, medium_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_baseline = Ivtot)

medium_mda_comp <- left_join(medium_baseline_Ivtot, medium_mda_Ivtot)

medium_mda_eff <- medium_mda_comp %>%
  summarise(diff_inf_vec = ((sum(Ivtot_baseline) - sum(Ivtot_int))/ sum(Ivtot_baseline))*100)

#compared to odin model
high_eir_odin_baseline_medium <- odin_baseline %>%
  filter(ref == 2 & between(t, medium_t_start, medium_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_baseline = Ivtot)

high_eir_odin_medium <- odin_mda_medium %>%
  filter(ref ==2 & between(t, medium_t_start, medium_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_int = Ivtot)

medium_mda_comp_odin <- left_join(high_eir_odin_baseline_medium, high_eir_odin_medium)

medium_mda_odin_eff <- medium_mda_comp_odin %>%
  summarise(diff_inf_vec = ((sum(Ivtot_baseline) - sum(Ivtot_int))/ sum(Ivtot_baseline))*100)

#eff diff in prevalence
medium_mda_prev <- msim_medium_bestfit %>%
  filter(t == medium_t_end+30&  endec_mu == 0.14 & wane == 0.01 ) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_int = slide_prev0to5)

baseline_prev_medium <- high_eir_msim_baseline %>%
  filter(t == medium_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_baseline = slide_prev0to5)

medium_mda_comp <- left_join(baseline_prev_medium, medium_mda_prev)
medium_mda_eff_prev <- medium_mda_comp %>%
  summarise(diff_prev = ((sum(prev_baseline) - sum(prev_int))/ sum(prev_baseline))*100)

#compared to odin model
high_eir_odin_prev <- odin_baseline %>%
  filter(ref == 2 & t == medium_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_baseline = slide_prev0to5)

high_eir_odin_medium_prev <- odin_mda_medium %>%
  filter(ref ==2 & t == medium_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_int = slide_prev0to5)

medium_mda_comp_odin <- left_join(high_eir_odin_prev, high_eir_odin_medium_prev)

medium_mda_odin_eff_prev <- medium_mda_comp_odin %>%
  summarise(diff_prev = ((sum(prev_baseline) - sum(prev_int))/ sum(prev_baseline))*100)

#eff diff inc u5
medium_mda_inc <- msim_medium_bestfit %>%
  filter(between(t, medium_t_start, medium_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_int = clin_inc0to5)

baseline_inc <- high_eir_msim_baseline %>%
  filter(between(t, medium_t_start, medium_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_baseline = clin_inc0to5)

medium_mda_comp <- left_join(baseline_inc, medium_mda_inc)
medium_mda_eff_inc <- medium_mda_comp %>%
  summarise(diff_inc = ((sum(inc_baseline) - sum(inc_int))/ sum(inc_baseline))*100)

#compared to odin model
high_eir_odin_inc <- odin_baseline %>%
  filter(ref == 2 & between(t, medium_t_start, medium_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_baseline = clin_inc0to5)

high_eir_odin_medium_inc <- odin_mda_medium %>%
  filter(ref ==2 & between(t, medium_t_start, medium_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_int = clin_inc0to5)

medium_mda_comp_odin <- left_join(high_eir_odin_inc, high_eir_odin_medium_inc)

medium_mda_odin_eff_inc <- medium_mda_comp_odin %>%
  summarise(diff_inc = ((sum(inc_baseline) - sum(inc_int))/ sum(inc_baseline))*100)

#repeat for early
early_mda_Ivtot <- msim_early_bestfit %>%
  filter(between(t, early_t_start, early_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_int = Ivtot)

early_baseline_Ivtot <-high_eir_msim_baseline %>%
  filter(between(t, early_t_start, early_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_baseline = Ivtot)

early_mda_comp <- left_join(early_baseline_Ivtot, early_mda_Ivtot)

early_mda_eff <- early_mda_comp %>%
  summarise(diff_inf_vec = ((sum(Ivtot_baseline) - sum(Ivtot_int))/ sum(Ivtot_baseline))*100)

#compared to odin model
high_eir_odin_baseline_early <- odin_baseline %>%
  filter(ref == 2 & between(t, early_t_start, early_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_baseline = Ivtot)

high_eir_odin_early <- odin_mda_early %>%
  filter(ref ==2 & between(t, early_t_start, early_t_end)) %>%
  select(t, Ivtot) %>%
  rename(Ivtot_int = Ivtot)

early_mda_comp_odin <- left_join(high_eir_odin_baseline_early, high_eir_odin_early)

early_mda_odin_eff <- early_mda_comp_odin %>%
  summarise(diff_inf_vec = ((sum(Ivtot_baseline) - sum(Ivtot_int))/ sum(Ivtot_baseline))*100)

#eff diff in prevalence
early_mda_prev <- msim_early_bestfit %>%
  filter(t == early_t_end+30&  endec_mu == 0.14 & wane == 0.01 ) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_int = slide_prev0to5)

baseline_prev_early <- high_eir_msim_baseline %>%
  filter(t == early_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_baseline = slide_prev0to5)

early_mda_comp <- left_join(baseline_prev_early, early_mda_prev)
early_mda_eff_prev <- early_mda_comp %>%
  summarise(diff_prev = ((sum(prev_baseline) - sum(prev_int))/ sum(prev_baseline))*100)

#compared to odin model
high_eir_odin_prev_early <- odin_baseline %>%
  filter(ref == 2 & t == early_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_baseline = slide_prev0to5)

high_eir_odin_early_prev <- odin_mda_early %>%
  filter(ref ==2 & t == early_t_end+30) %>%
  select(t, slide_prev0to5) %>%
  rename(prev_int = slide_prev0to5)

early_mda_comp_odin <- left_join(high_eir_odin_prev_early, high_eir_odin_early_prev)

early_mda_odin_eff_prev <- early_mda_comp_odin %>%
  summarise(diff_prev = ((sum(prev_baseline) - sum(prev_int))/ sum(prev_baseline))*100)

#eff diff inc u5
early_mda_inc <- msim_early_bestfit %>%
  filter(between(t, early_t_start, early_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_int = clin_inc0to5)

baseline_inc_early <- high_eir_msim_baseline %>%
  filter(between(t, early_t_start, early_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_baseline = clin_inc0to5)

early_mda_comp <- left_join(baseline_inc_early, early_mda_inc)
early_mda_eff_inc <- early_mda_comp %>%
  summarise(diff_inc = ((sum(inc_baseline) - sum(inc_int))/ sum(inc_baseline))*100)

#compared to odin model
high_eir_odin_inc <- odin_baseline %>%
  filter(ref == 2 & between(t, early_t_start, early_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_baseline = clin_inc0to5)

high_eir_odin_early_inc <- odin_mda_early %>%
  filter(ref ==2 & between(t, early_t_start, early_t_start+180)) %>%
  select(t, clin_inc0to5) %>%
  rename(inc_int = clin_inc0to5)

early_mda_comp_odin <- left_join(high_eir_odin_inc, high_eir_odin_early_inc)

early_mda_odin_eff_inc <- early_mda_comp_odin %>%
  summarise(diff_inc = ((sum(inc_baseline) - sum(inc_int))/ sum(inc_baseline))*100)

timing_mda <- c("late", "medium", "early")
endpoint <- c("Ivtot", "slide_prev0to5", "clin_inc0to5")

late_Ivtot <- data.frame(late_mda_eff, late_mda_odin_eff)
names(late_Ivtot) <- c("msim_mda", "odin_mda")
late_Ivtot <- late_Ivtot %>%
  mutate(timing_mda = timing_mda[1],
         endpoint = endpoint[1])

late_prev <- data.frame(late_mda_eff_prev, late_mda_odin_eff_prev)
names(late_prev) <- c("msim_mda", "odin_mda")
late_prev <- late_prev %>%
  mutate(timing_mda = timing_mda[1],
         endpoint = endpoint[2])

late_inc <- data.frame(late_mda_eff_inc, late_mda_odin_eff_inc)
names(late_inc) <- c("msim_mda", "odin_mda")
late_inc <- late_inc %>%
  mutate(timing_mda = timing_mda[1],
         endpoint = endpoint[3])

medium_Ivtot <- data.frame(medium_mda_eff, medium_mda_odin_eff)
names(medium_Ivtot) <- c("msim_mda", "odin_mda")
medium_Ivtot <- medium_Ivtot %>%
  mutate(timing_mda = timing_mda[2],
         endpoint = endpoint[1])

medium_prev <- data.frame(medium_mda_eff_prev, medium_mda_odin_eff_prev)
names(medium_prev) <- c("msim_mda", "odin_mda")
medium_prev <- medium_prev %>%
  mutate(timing_mda = timing_mda[2],
         endpoint = endpoint[2])

medium_inc <- data.frame(medium_mda_eff_inc, medium_mda_odin_eff_inc)
names(medium_inc) <- c("msim_mda", "odin_mda")
medium_inc <- medium_inc %>%
  mutate(timing_mda = timing_mda[2],
         endpoint = endpoint[3])

early_Ivtot <- data.frame(early_mda_eff, early_mda_odin_eff)
names(early_Ivtot) <- c("msim_mda", "odin_mda")
early_Ivtot <- early_Ivtot %>%
  mutate(timing_mda = timing_mda[3],
         endpoint = endpoint[1])

early_prev <- data.frame(early_mda_eff_prev,early_mda_odin_eff_prev)
names(early_prev) <- c("msim_mda", "odin_mda")
early_prev <- early_prev %>%
  mutate(timing_mda = timing_mda[3],
         endpoint = endpoint[2])

early_inc <- data.frame(early_mda_eff_inc, early_mda_odin_eff_inc)
names(early_inc) <- c("msim_mda", "odin_mda")
early_inc <- early_inc %>%
  mutate(timing_mda = timing_mda[3],
         endpoint = endpoint[3])

mda_timing_df <- do.call("rbind", list(late_Ivtot, late_prev, late_inc,
                                       medium_Ivtot, medium_prev, medium_inc,
                                       early_Ivtot, early_prev, early_inc))

mda_timing_df_long_high_EIR <- mda_timing_df %>%
  pivot_longer(cols = c("msim_mda", "odin_mda"), names_to = "model", values_to = "eff") %>%
  mutate(EIR = "high")

mda_timing_df_long_low_EIR <- readRDS(file = "analysis/exploring_interactions/MIM_poster/chapter_plots/mda_timing_df_long_low_EIR.rds")

mda_timing_df_long <- rbind(mda_timing_df_long_high_EIR, mda_timing_df_long_low_EIR)


#plot comparing efficacy of the different models
mda_timing_df_long$model <- factor(mda_timing_df_long$model,
                                   levels = c("odin_mda", "msim_mda"))


sensitivity_analysis_model_fits <- ggplot(mda_timing_df_long,
                                          aes(x = factor(timing_mda, levels = c("early", "medium", "late")),
                                              y = eff, fill = as.factor(model))) +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  facet_grid(EIR ~ endpoint,
             labeller = labeller(
               endpoint = c(
                 slide_prev0to5 = "Efficacy (%) on \n slide prevalence in under 5-year-olds",
                 clin_inc0to5  = "Cases averted (%) \n in under 5-year-olds",
                 Ivtot         = "Reduction (%) in \n infectious vectors"
               ),
               EIR = c(
                 low  = "Low transmission setting",
                 high = "High transmission setting"
               )
             )) +
  theme_bw(base_size = 14) +
  labs(x = "Timing of MDA in relation to ITN campaign", y = "Efficacy (%)") +
  theme(legend.position = "bottom") +
  ylim(0, 100) +
  scale_fill_manual(name = "Model",
                    values = c('#1b9e77', '#d95f02'),
                    labels = c("Model A", "Model C"))

ggsave(sensitivity_analysis_model_fits, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/odin_sensitivity_plots_mda_timing_eir.pdf")

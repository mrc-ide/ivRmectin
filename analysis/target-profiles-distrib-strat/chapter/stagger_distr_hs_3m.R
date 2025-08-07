#model staggered distribution by altering hazard ratios.
#devtools::load_all()
#we assume equal coverage at each distribution

#modify the hazard ratio curve
require(tidyverse)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")


stag_10 <- 11
times_10_d <- seq(1, stag_10, length.out = 3)

ivm_haz <- ivm_haz %>%
  select(-IVM_400_1_HS) %>%
  filter(between(Day, 1, 23))

ivm_haz_df <- as.data.frame(ivm_haz)

shift_days_10_MDA_3 <- times_10_d[3] - times_10_d[1]

max_day <- max(ivm_haz_df$Day)
new_days <- (max_day + 1):90 # whole thing will be completed within 90 days

extra_rows_10d <- data.frame(Day = new_days, IVM_300_3_HS = 1)

df_extended_10d <- rbind(ivm_haz_df, extra_rows_10d)

group1_HR <- c(df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7))

df_extended_10d <- df_extended_10d %>%
  mutate(group1 = group1_HR)

#create new shifted column - group 3 - and fill with HR 1 first, no impact if not started yet

df_extended_10d$group3 <- 1 # HR1 - no impact

df_extended_10d$group3[(shift_days_10_MDA_3 + 1):nrow(df_extended_10d)] <- df_extended_10d$group1[1:(nrow(df_extended_10d) - shift_days_10_MDA_3)]

shift_days_10_MDA_2 <- times_10_d[2] - times_10_d[1]

df_extended_10d$group2 <- 1
df_extended_10d$group2[(shift_days_10_MDA_2 + 1):nrow(df_extended_10d)] <- df_extended_10d$group1[1:(nrow(df_extended_10d) - shift_days_10_MDA_2)]

df_extended_10d <- df_extended_10d %>%
  select(Day, IVM_300_3_HS, group1, group2, group3)

#df_extended_10d_new_HR <- df_extended_10d %>%
#  #mutate(across(everything(), ~ replace_na(.x, 1))) %>% #no killing effect of others as MDA not started yet
#  rowwise() %>%
#  mutate(HR_use = mean(c(group1, group2, group3))) %>%
#  mutate(stagger = "10d")

df_extended_10d_new_HR <- df_extended_10d %>%
  rowwise() %>%
  mutate(HR_use_above1 = mean(c_across(group1:group3)[c_across(group1:group3) > 1], na.rm = TRUE), #rowwise mean, if greater than 1
         prop_pop_cov = sum(c_across(group1:group3) > 1)/3,
         stagger = "10d") %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  ungroup()

#for the 20d strategy

stag_20 <- 21
times_20_d <- seq(1, stag_20, length.out = 3)

shift_days_20_MDA_3 <- times_20_d[3] - times_20_d[1]

max_day <- max(ivm_haz_df$Day)
new_days <- (max_day + 1):123

extra_rows_20d <- data.frame(Day = new_days, IVM_300_3_HS = 1)

df_extended_20d <- rbind(ivm_haz_df, extra_rows_20d)

to_add <- nrow(df_extended_20d) - 90

group1_HR <- c(df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               rep(1, to_add))

length(group1_HR)

df_extended_20d <- df_extended_20d %>%
  mutate(group1 = group1_HR)

#create new shifted column - group 3 - and fill with HR 1 first, no impact

df_extended_20d$group3 <- 1 #HR 1 - no impact

df_extended_20d$group3[(shift_days_20_MDA_3 + 1):nrow(df_extended_20d)] <- df_extended_20d$group1[1:(nrow(df_extended_20d) - shift_days_20_MDA_3)]

shift_days_20_MDA_2 <- times_20_d[2] - times_20_d[1]

df_extended_20d$group2 <- 1
df_extended_20d$group2[(shift_days_20_MDA_2 + 1):nrow(df_extended_20d)] <- df_extended_20d$group1[1:(nrow(df_extended_20d) - shift_days_20_MDA_2)]

df_extended_20d <- df_extended_20d %>%
  select(Day, IVM_300_3_HS, group1, group2, group3)

df_extended_20d <- df_extended_20d[1:113,]

df_extended_20d_new_HR <- df_extended_20d %>%
  rowwise() %>%
  mutate(HR_use_above1 = mean(c_across(group1:group3)[c_across(group1:group3) > 1], na.rm = TRUE), #rowwise mean, if greater than 1
         prop_pop_cov = sum(c_across(group1:group3) > 1)/3,
         stagger = "20d") %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  ungroup()


HR_all <- rbind(df_extended_10d_new_HR, df_extended_20d_new_HR)

head(HR_all)

ggplot(HR_all, aes(x = Day, y = group1))+
  geom_line()+
  geom_line(aes(x = Day, y = group2,col = "red"))+
  geom_line(aes(x = Day, y = group3, col = "green"))+
  geom_line(aes(x = Day, y = HR_use_above1), col = "orange")+
  facet_wrap(vars(stagger))+
  theme_bw()

HR_all_long <- HR_all %>%
  pivot_longer(cols = c(group1, group2, group3, HR_use_above1), names_to = "group",
               values_to = "HR") %>%
  mutate(group_lab = case_when(group == "HR_use_above1" ~ "average",
                               TRUE ~ "group-level"),
         group_combo = paste0(group_lab))

HR_plot <- ggplot(HR_all_long, aes(x = Day, y = HR, group = group))+
  geom_line(aes(col = as.factor(group), lty = as.factor(group_lab)),size = 1.1,
            inherit.aes = TRUE)+
  geom_point(aes(col = as.factor(group)))+
  facet_wrap(vars(stagger), labeller = label_both)+
  theme_bw()+
  geom_hline(aes(yintercept = 1), lty = "dashed")+
  ylim(1, 10)+
  scale_colour_manual(values = c(
    "group1" = "#1b9e77",
    "group2" = "#d95f02",
    "group3" = "#7570b3",
    "HR_use_above1" = "black"), labels= c("Group 1", "Group 2", "Group 3", "Group average"),
    name = "Hazard ratio group"
  )+
  labs(linetype = "Average or group-level", y = "Hazard ratio")+
  ggtitle("Staggered distribution")


time_cov_plot <- ggplot(HR_all, aes(x = Day, y = prop_pop_cov))+
  geom_line()+
  facet_wrap(vars(stagger), labeller = label_both)+
  ylab("Proportion of covered group with HR > 1")+
  theme_bw()+
  ylim(0,1)

write_rds(HR_all_long, file = "analysis/target-profiles-distrib-strat/chapter/output/HR_staggered_long.rds")
write_rds(HR_all, file = "analysis/target-profiles-distrib-strat/chapter/output/HR_staggered.rds")

HR_all_in_one_plot <- ggplot(ivm_haz, aes(x = Day, y = IVM_300_3_HS))+
  geom_line()+
  geom_point()+
  theme_bw()+
  ylim(0,10)+
  ylab("Hazard ratio")+
  ggtitle("All in one")



staggered_HR <- cowplot::plot_grid(HR_plot, time_cov_plot, align = "v",
                   nrow = 2)

HR_figure <- cowplot::plot_grid(HR_all_in_one_plot, staggered_HR)



df_all <-df_extended_10d_new_HR %>%
  select(Day, group1) %>%
  rename(HR_use_above1 = group1) %>% #giving this name so easy to bind with others
  mutate(HR_use_above1 = case_when(HR_use_above1 == 1 ~ round(1/1),
                                     TRUE ~ HR_use_above1)) %>%
  mutate(stagger = "all_in_one") %>%
  rowwise() %>%
  mutate(
         prop_pop_cov = sum(c_across(HR_use_above1) > 1)/1)



HR_plot_all_1 <- ggplot(df_all, aes(x = Day, y = HR_use_above1))+
  geom_line()+
  theme_minimal()+
  ylab("Hazard Ratio")

time_cov_all_1 <- ggplot(df_all, aes(x = Day, y = prop_pop_cov))+
  geom_line()+
  theme_minimal()+
  ylab("prop_pop_cov")

cowplot::plot_grid(HR_plot_all_1, time_cov_all_1, align = "v",
                   nrow = 2)

distr_pals <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')

#HR_all_plot <- HR_all %>%
#  select(Day, group1, HR_use_above1, stagger)
#
#HR_all_in_one <- HR_all_plot %>%
#  ungroup() %>%
#  select(Day, HR_use_above1 = group1) %>%
#  mutate(stagger = "all-in-one")
#
#HR_combined <- bind_rows(HR_all_plot, HR_all_in_one) %>%
#  select(-group1) %>%
#  mutate(stagger = factor(stagger, levels = c("all-in-one", "10d", "20d")))
#
#
#
#HR_plot <- ggplot(HR_combined, aes(x = Day, y = HR_use_above1, col = as.factor(stagger)))+
#  geom_line(size = 1.1)+
#  theme_bw()+
#  scale_colour_manual(values = c(distr_pals[3], distr_pals[1], distr_pals[2]), name = "Time to complete 1 round of MDA",
#                      labels = c("1 day", "10 days", "20 days"))+
#  ylim(0, 10)+
#  ylab("Mean hazard ratio")+
#  xlim(1, 113)+
#  theme(legend.position = c(0.8, 0.8))+
#  geom_hline(aes(yintercept = 1), lty = "dashed")
#ggsave(HR_plot, file = "analysis/target-profiles-distrib-strat/chapter/plots/HR_inputs.svg")
#ggsave(HR_plot, file = "analysis/target-profiles-distrib-strat/chapter/plots/HR_inputs.pdf")

#see clear peaks for the 10d, the period of overlap within each MDA is big but between months not much.
#plot this.
#for 30d, month 1 and month 2, or month 2 and 3 overlap because of the staggering across groups is bigger.

#ggplot(df_extended_10d_new_HR, aes(x = Day, y = group1), col = "black")+
#  geom_line()+
#  geom_line(aes(x = Day, y = group2), col = "red")+
#  geom_line(aes(x = Day, y = group3), col = "blue")
#
#ggplot(df_extended_20d_new_HR, aes(x = Day, y = group1), col = "black")+
#  geom_line()+
#  geom_line(aes(x = Day, y = group2), col = "red")+
#  geom_line(aes(x = Day, y = group3), col = "blue")


#we now pass this into the model
devtools::load_all()
time_period <- 365*15 #long run to get to eqm
mda_int <- 30
#ivm_cov_in = c(0.1, 0.9)

# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

#itn_on <- 100 #introduce nets 100 days into simulation

#modelling with extreme density dependence

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

#distr plan
start <- (365*5)+200 #time for starting distr, help get to reach eqm, and to distr at right time when seasonality introduced

#Q0_in <- c(0.21, 0.71, 0.92, 0.94)
Q0_in <- 0.92
init_EIR_in <- c(2, 100)
#init_EIR_in <- 100
#target_cov <- c(0.6, 0.7, 0.8)
target_cov <- 0.7

df_var_all <- expand.grid(Q0 = Q0_in, ivm_cov = target_cov, init_EIR = init_EIR_in)
#df <- df[1,]
my_list_all <- list()
for (i in seq_len(nrow(df_var_all))){
  my_list_all[[i]] <- as.numeric(df_var_all[i,])
}



ivm_parms_10d_stag <- ivRmectin::ivm_fun_stag_cov(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  prop_human_HR_threshold = df_extended_10d_new_HR$prop_pop_cov,
  time_period = time_period,
  hazard_profile = df_extended_10d_new_HR$HR_use_above1,
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

range(ivm_parms_10d_stag$prop_human_HR_threshold)
plot(ivm_parms_10d_stag$prop_human_HR_threshold)


mod_10d_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_10d_stag$ttt,
    eff_len = ivm_parms_10d_stag$eff_len,
    haz = ivm_parms_10d_stag$haz,
    prop_human_HR_threshold = ivm_parms_10d_stag$prop_human_HR_threshold,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_10d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_10d_stag$ivm_max_age,
    IVRM_start = ivm_parms_10d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_10d_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_10d_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0,  ivm_cov, slide_prev0to5, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "10d-stagger"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_10d_stag <- my_sim_mod_10d_stag()

ggplot(df_mod_10d_stag, aes(x = t, y = ivm_cov))+
  geom_line()

max(df_mod_10d_stag$ivm_cov) #only short period of time when achieve that coverage


#summarise the coverage outputs
df_mod_10d_stag %>%
  filter(between(t, start, start + 23)) %>%
  summarise(mean_cov = mean(ivm_cov),
            median_cov = median(ivm_cov),
            min_cov = min(ivm_cov),
            max_cov = max(ivm_cov))

#then 20d stagger
ivm_parms_20d_stag <- ivRmectin::ivm_fun_stag_cov(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  prop_human_HR_threshold = df_extended_20d_new_HR$prop_pop_cov,
  time_period = time_period,
  hazard_profile = df_extended_20d_new_HR$HR_use_above1,
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

mod_20d_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_20d_stag$ttt,
    prop_human_HR_threshold = ivm_parms_20d_stag$prop_human_HR_threshold,
    eff_len = ivm_parms_20d_stag$eff_len,
    haz = ivm_parms_20d_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_20d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_20d_stag$ivm_max_age,
    IVRM_start = ivm_parms_20d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_20d_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_20d_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "20d-stagger"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_20d_stag <- my_sim_mod_20d_stag()

#all-in-one
ivm_parms_all <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = c(start, start + 30, start + 60),
  #IVM_start_times = start,
  time_period = time_period,
  hazard_profile = df_all$HR_use_above1[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

mod_all <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_all <- function(){
  mod_out_list <- lapply(my_list_all, mod_all)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, Svtot, Evtot, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "all-in-one"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_all <- my_sim_mod_all()

#checking against the all-in-one in the staggered model
ivm_parms_all_stag <- ivRmectin::ivm_fun_stag_cov(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  time_period = time_period,
  hazard_profile = df_all$HR_use_above1,
  prop_human_HR_threshold = df_all$prop_pop_cov, #still passing into ivm compartment but experience HR of 1
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)


df_all <- df_all %>%
  mutate(t = Day+start-1)

ivm_parms_all_stag$eff_len

ggplot(df_all, aes(x = Day, y = HR_use_above1))+
  geom_line()+
  geom_point()

mod_all_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_all_stag$ttt,
    eff_len = ivm_parms_all_stag$eff_len,
    prop_human_HR_threshold = ivm_parms_all_stag$prop_human_HR_threshold,
    haz = ivm_parms_all_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_all_stag$ivm_min_age,
    ivm_max_age = ivm_parms_all_stag$ivm_max_age,
    IVRM_start = ivm_parms_all_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_all_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_all_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5,Svtot, Evtot, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "all-in-one-stag"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_all_stag <- my_sim_mod_all_stag()

#compare stag and normal code, for just one distribution (i.e. not actually any staggering.)
df_mod_compare <- rbind(df_mod_all, df_mod_all_stag)

df_mod_compare <- df_mod_compare %>%
  mutate(init_EIR = case_when(ref == 1 ~ 2,
                             TRUE ~ 100)) %>%
  filter(init_EIR == 100)

ggplot(df_mod_compare, aes(x = t, y = 1-ivm_cov, col = as.factor(model_type)))+
  geom_line()

ggplot(df_mod_compare, aes(x = t, y = ivm_cov, col = as.factor(model_type)))+
  geom_line()+
  xlim(2000, 2200)+
  ylab("Ivermectin coverage")

ggplot(df_mod_compare, aes(x = t, y = (ivm_cov + (1-ivm_cov)), col = as.factor(model_type)))+
  geom_line()

ggplot(df_mod_compare, aes(x = t, y = mv, col = as.factor(model_type)))+
  geom_line()+
  theme_minimal()+
  xlim(2000, 2200)+
  ylim(20,45)

ggplot(df_mod_compare, aes(x = t, y = Ivtot, col = as.factor(model_type)))+
  geom_line()+
  theme_minimal()+
  xlim(2000, 2200)

ggplot(df_mod_compare, aes(x = t, y = Svtot, col = as.factor(model_type)))+
  geom_line()+
  theme_minimal()+
  xlim(2000, 2200)

ggplot(df_mod_compare, aes(x = t, y = Evtot, col = as.factor(model_type)))+
  geom_line()+
  theme_minimal()+
  xlim(2000, 2200)

ggplot(df_mod_compare, aes(x = t, y = EIRout, col = as.factor(model_type)))+
  geom_line()+
  theme_minimal()+
  xlim(2000, 2200)

ggplot(df_mod_compare, aes(x = t, y = slide_prev0to5*100, col = as.factor(model_type)))+
  geom_line()+
  theme_minimal()+
  xlim(2000, 2500)+
  ylim(0,100)+
  ylab("Slide prevalence (%) under 5s")

haz_check_all <- ivm_parms_all$haz
haz_check_stag <- ivm_parms_all_stag$haz


haz_check_all == haz_check_stag[1:23]

haz_check_all == haz_check_stag[31:53]

haz_check_all == haz_check_stag[61:83]


df_mod_compare %>%
  filter(between(t,start, start+90))%>%
  group_by(model_type)%>%
  summarise(tot_mosq = sum(mv),
            mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5),
            tot_Iv = sum(Ivtot))


ggplot(df_mod_compare, aes(x = t, y = clin_inc0to5*1000, col = as.factor(model_type)))+
  geom_line()+
  ylim(0, 6)+
  xlim(start, start + 90)

ggplot(df_mod_compare, aes(x = t, y = slide_prev0to5, col = as.factor(model_type)))+
  geom_line()+
  ylim(0, 1)+
  xlim(2000, 2500)

ggplot(df_mod_compare, aes(x = t, y = mv, col = as.factor(model_type)))+
  geom_line()+
  xlim(start, start+100)+
  geom_point()
#baseline scenario

mod_baseline <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = 0, #make coverage 0
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_baseline <- function(){
  mod_out_list <- lapply(my_list_all, mod_baseline)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "baseline"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_baseline <- my_sim_mod_baseline()

df_mod_distr <- do.call("rbind", list(df_mod_10d_stag, df_mod_20d_stag, df_mod_all, df_mod_baseline))

write_rds(df_mod_distr, file = "analysis/target-profiles-distrib-strat/chapter/output/df_distr_HR_3m_HS.rds")

#with seasonality

#10d Stagger

mod_10d_stag_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_10d_stag$ttt,
    eff_len = ivm_parms_10d_stag$eff_len,
    haz = ivm_parms_10d_stag$haz,
    prop_human_HR_threshold = ivm_parms_10d_stag$prop_human_HR_threshold,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_10d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_10d_stag$ivm_max_age,
    IVRM_start = ivm_parms_10d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_10d_stag_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_10d_stag_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "10d-stagger-Sen"))}, simplify = F))
  return(mod_df)
}

df_mod_10d_stag_Sen <- my_sim_mod_10d_stag_Sen()

#over 30 days

mod_20d_stag_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_20d_stag$ttt,
    eff_len = ivm_parms_20d_stag$eff_len,
    haz = ivm_parms_20d_stag$haz,
    prop_human_HR_threshold = ivm_parms_20d_stag$prop_human_HR_threshold,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_20d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_20d_stag$ivm_max_age,
    IVRM_start = ivm_parms_20d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_20d_stag_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_20d_stag_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "20d-stagger-Sen"))}, simplify = F))
  return(mod_df)
}

df_mod_20d_stag_Sen <- my_sim_mod_20d_stag_Sen()

#all in one, with seasonality
mod_all_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_all_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_all_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "all-in-one-Sen"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_all_Sen <- my_sim_mod_all_Sen()

#baseline
mod_baseline_Sen <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = 0, #make coverage 0
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_baseline_Sen <- function(){
  mod_out_list <- lapply(my_list_all, mod_baseline_Sen)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "baseline-Sen"))}, simplify = F))
  return(mod_df)
}

df_mod_baseline_Sen <- my_sim_mod_baseline_Sen()

df_mod_distr_Sen <- do.call("rbind", list(df_mod_10d_stag_Sen, df_mod_20d_stag_Sen, df_mod_all_Sen, df_mod_baseline_Sen))

write_rds(df_mod_distr_Sen, file = "analysis/target-profiles-distrib-strat/chapter/output/df_distr_HR_Sen_3m_HS.rds")

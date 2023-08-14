#Do LLINs influence the excess mortality rate due to ivermectin#?
#going to use this analysis to see if we need to change the way in which we model excess mortality due to ivm in malariasimulation

#odin_endec_mu_h_model.R has the mu_vi[i] = haz[i]*(mu+mu_h), where haz is always 1

#odin_endec_mu_h_constant_ivm_uptake.R #this model has a constant parameter for the uptake of ivermectin

require(tidyverse)


#script for running endec_mosq_model
# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
# Create a vector of age categories for the model
#init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 100 #low - 2, moderate - 15, high - 120

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years
time_period <- 365*10
IVM_begin <- (365*8)+180
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)
ivm_cov = 0.8

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}


#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

eff_len <- 23
ivm_on <- IVM_start[1]
ivm_off <- IVM_start[3]+eff_len

ivm_haz2 = rep(1, time_period) #not using hazards, so just set hazards to 1
ivm_parms3 <- ivm_fun(IVM_start_times = IVM_start,# time endectocide delivery occurs
                      #IVM_start = 180,
                      time_period = time_period,         # time period for the model to run over
                      hazard_profile = ivm_haz2[1:23], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                      #hazard_profile = ivm_haz2[1:730],
                      ivm_coverage = ivm_cov, # proportion of population receiving the endectocide
                      ivm_min_age = 5, # youngest age group receiving endectocide
                      ivm_max_age = 90) # oldest age group receiving endectocide

#from Ellie's work
pyr_only_dn0 <- 0.341
pyr_only_rn0 <- 0.637

#run Hannah's model, output the avhc and mosquito density over time
wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  #num_int = 1,
                                  num_int = 1,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model
res0 <- runfun(wh0)
df_0 <- as.data.frame(res0)

avhc_constant <- unique(df_0$avhc) #0.3066667
#avhc_constant <- 0.3066667

#plot mosquito density over time
ggplot(df_0, aes(x = t, y = mv))+
  geom_line()

#from d180 to d240+23
df_0_ivm <- df_0 %>%
  filter(between(t, ivm_on, ivm_off))

#mv during the ivm distrib time
#then select mv during the ivermectin distribution time




#then run the simplified version for different values of mu_h and see which values gives best fit to Hannah's mosquito densities
create_mu_h_loop <- function(mu_h_in) {
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_endec_mu_h_model.R",
    num_int = 1,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms3$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms3$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms3$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_parms3$ivm_cov_par, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms3$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms3$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms3$IVRM_start,
    mu_h = mu_h_in
  )
  return(output)
}

mu_h_vector <- seq(0, 1, 0.01)

#generate parameter set
out_lapply_list <- lapply(mu_h_vector, create_mu_h_loop)

#the run the model
res_out_list <- lapply(out_lapply_list, runfun)

#save as a df

df_1 <- as.data.frame(res_out_list)
#write.csv(df_1, file = "data/llin_ivm_muh/df_1.csv")




res_ivm_distrib <- df_1 %>%
  filter(between(t, ivm_on, ivm_off))


#go through the list form output and save t, mv, avhc and mu_h
df_1a <- do.call(rbind,
                  sapply(1:length(mu_h_vector), function(x){
                    as.data.frame(res_out_list[[x]]) %>% #go through all the values in the list and save
                      filter(between(t, ivm_on, ivm_off)) %>%
                      select(t, mv, avhc, mu_h) %>%
                      mutate(ref=x)
                  }, simplify = F))

out_list <- split(df_1a, f = df_1a$ref)

#compute sum squared error

error <- numeric()
for (i in 1:length(mu_h_vector)){
  error <- c(error, sum((df_0_ivm $mv - out_list[[i]]$mv)^2))
}
error
index <- which.min(error) #index the smallest value
mu_h_fit <- mu_h_vector[index]

#mu_h_fit is 0.26

#put this back into the odin_endec_mu_h_model.R and check fit to Hannah's model

#ivm only, no nets
wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_endec_mu_h_model.R",
                                  #num_int = 1,
                                  num_int = 1,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms3$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  mu_h = 0.26,
                                  IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model with mu_h of 0.26
res2 <- runfun(wh2)
df_2 <- as.data.frame(res2)

#NC model, nets only, no ivm
wh2a <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_endec_mu_h_model.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms3$haz,
                                  ivm_cov_par = 0,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  mu_h = 0.26,
                                  IVRM_start = ivm_parms1$IVRM_start,
                                  d_ITN0 = pyr_only_dn0,
                                  r_ITN0 = pyr_only_rn0)

#run Hannah's model with mu_h of 0.26
res2a <- runfun(wh2a)
df_2a <- as.data.frame(res2a)

#can see that mv broadly matches
mv_HS_NC <- ggplot(df_0, aes(x = t, y = mv))+ #HS
  geom_line()+
  geom_line(data = df_2, aes(x = t, y = mv, col = "red"), show.legend = FALSE)+ #NC
  theme_minimal()+
  ggtitle("mv w HS model (black) and NC model (red)")+
  xlim(2920, 3650)
#ggsave(mv_HS_NC, file = "plots/llin_ivm_muh/mv_HS_NC.png")

#rbind them
df_0_bind <- df_0 %>%
  mutate(model = "HS") %>%
  select(t, mv, model)

df_2_bind <- df_2 %>%
  mutate(model = "NC") %>%
  select(t, mv, model)

mv_fit_data <- rbind(df_0_bind, df_2_bind)
write.csv(mv_fit_data, file = "data/llin_ivm_muh/mv_fit_data.csv", row.names = FALSE)

breaks_plot <- seq(2920, 3285, by = 30)

mv_plot_esa <- ggplot(mv_fit_data, aes(x = t, y = mv, col = as.factor(model)))+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "palegreen", alpha =0.25)+
  geom_line()+
  #xlim(2920, 3650)+
  ylab("Mosquito Density")+
  theme_minimal()+
  ylim(0, 50)+
  scale_x_continuous(limits = c(2920, 3285), breaks = breaks_plot, minor_breaks = breaks_plot, name = "Time (days)")+
  scale_colour_manual(values = c("black", "blue"), labels = c("Slater et al (2020) model: \n ivermectin mortality rate modelled with daily hazard ratios", "Modified Slater et al (2020) model: \n fixed additional mortality rate during ivermectin distribution period"),
                      name = "Modelling method")+
  theme(legend.position = c(.25, .95))+
  theme(text = element_text(size = 14))+
  #geom_vline(xintercept = c(ivm_on, ivm_off), linetype = "dashed")+
  annotate("segment", x = IVM_start[1], xend = IVM_start[1], y = 45, yend = 42, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[2], xend = IVM_start[2], y = 45, yend = 42, colour = "red", arrow = arrow())+
  annotate("segment", x = IVM_start[3], xend = IVM_start[3], y = 45, yend = 42, colour = "red", arrow = arrow())
  #geom_rect(aes(xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf), fill = "green", alpha = 0.05)
#ggsave(mv_plot_esa, file = "plots/llin_ivm_muh/mv_plot_esa.svg)


#then run Hannah's models with LLINs (cov = 60%) and output the avhc
wh3 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start,
                                  d_ITN0 = pyr_only_dn0,
                                  r_ITN0 = pyr_only_rn0)

#run Hannah's model
res3 <- runfun(wh3)
df_3 <- as.data.frame(res3)
unique(df_3$avhc)

#nets only, no ivm
wh3a <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = 0,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start,
                                  d_ITN0 = pyr_only_dn0,
                                  r_ITN0 = pyr_only_rn0)

#run Hannah's model
res3a <- runfun(wh3a)
df_3a <- as.data.frame(res3a)

#untreated nets and ivm. param's from Ettie's work
wh5 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start,
                                  r_ITN0 = 0.409,
                                  d_ITN0 = 0.059)

#run Hannah's model
res5 <- runfun(wh5)
df_5 <- as.data.frame(res5)


#then put into the model with constant ivm uptake, and we will sub in the mu_H value that we just estimated
wh4 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_endec_mu_h_constant_ivm_uptake.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms3$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  mu_h = 0.26,
                                  avhc = avhc_constant,
                                  IVRM_start = ivm_parms1$IVRM_start,
                                  d_ITN0 = pyr_only_dn0,
                                  r_ITN0 = pyr_only_rn0)

#run Hannah's model with mu_h of 0.26
res4 <- runfun(wh4)
df_4 <- as.data.frame(res4)


#NC model with untreated net: dn0 = 0.059, rn0 = 0.409, sn0 = 0.532
wh6 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_endec_mu_h_constant_ivm_uptake.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms3$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  mu_h = 0.26,
                                  avhc = avhc_constant,
                                  IVRM_start = ivm_parms1$IVRM_start,
                                  d_ITN0 = 0.059,
                                  r_ITN0 = 0.409)

#run Hannah's model with mu_h of 0.26
res6 <- runfun(wh6)
df_6 <- as.data.frame(res6)



#df_0 : Hannah's model, no nets, ivm
df_0_epi <- df_0 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead)
write.csv(df_0_epi, file = "data/llin_ivm_muh/df_0_epi.csv", row.names = FALSE)

#df_3: Hannah's model, nets, ivm
df_3_epi <- df_3 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm,mv_dead, mvx_dead)
write.csv(df_3_epi, file = "data/llin_ivm_muh/df_3_epi.csv", row.names = FALSE)

#df_3a: Hannah's model, nets, no ivm
df_3a_epi <- df_3a %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm,mv_dead, mvx_dead)
write.csv(df_3a_epi, file = "data/llin_ivm_muh/df_3a_epi.csv", row.names = FALSE)

#df_5: Hannah's model, untreated nets, ivm
df_5_epi <- df_5 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm,mv_dead, mvx_dead)
write.csv(df_5_epi, file = "data/llin_ivm_muh/df_5_epi.csv", row.names = FALSE)

#df2: Nilani's model, no nets, ivm
df_2_epi <- df_2 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead)
write.csv(df_2_epi, file = "data/llin_ivm_muh/df_2_epi.csv", row.names = FALSE)

#df2a: Nilani's model, nets only
df_2a_epi <- df_2a %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead)
write.csv(df_2a_epi, file = "data/llin_ivm_muh/df_2a_epi.csv", row.names = FALSE)


#df4: Nilani's model, nets, ivm
df_4_epi <- df_4 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead)
write.csv(df_4_epi, file = "data/llin_ivm_muh/df_4_epi.csv", row.names = FALSE)

#df6: NC model, untreated nets, ivm

df_6_epi <- df_6 %>%
  select(t, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead)
write.csv(df_6_epi, file = "data/llin_ivm_muh/df_6_epi.csv", row.names = FALSE)


#READ IN THE EPI FILES HERE

#HS models

#HS, IVM only
df_0_epi <- read.csv("data/llin_ivm_muh/df_0_epi.csv", header = TRUE)

#HS, IVM + NETS
df_3_epi <- read.csv("data/llin_ivm_muh/df_3_epi.csv", header = TRUE)


#HS, NETS only
df_3a_epi <- read.csv("data/llin_ivm_muh/df_3a_epi.csv", header = TRUE)

#HS, UTN and ivm
df_5_epi <- read.csv("data/llin_ivm_muh/df_5_epi.csv", header = TRUE)


#NC models

#NC, IVM only
df_2_epi <- read.csv("data/llin_ivm_muh/df_2_epi.csv", header = TRUE)


#NC, IVM + NETS
df_4_epi <- read.csv("data/llin_ivm_muh/df_4_epi.csv", header = TRUE)


#NC, NETS only
df_2a_epi <- read.csv("data/llin_ivm_muh/df_2a_epi.csv", header = TRUE)

#NC, IVM + UTN
df_6_epi <- read.csv("data/llin_ivm_muh/df_6_epi.csv", header = TRUE)


#plots to compare the methods
NC_HS_model_nets <- ggplot(df_3_epi, aes(x = t, y = EIR_tot))+ #HS
  geom_line()+
  geom_line(data = df_4_epi, aes(x = t, y = EIR_tot, col = "red"), show.legend = FALSE)+ #NC
  theme_minimal()+
  ggtitle("HS model (black), NC model (red), w nets")

NC_HS_model_no_nets <-ggplot(df_0_epi, aes(x = t, y = EIR_tot))+ #HS
  geom_line()+
  geom_line(data = df_2_epi, aes(x = t, y = EIR_tot, col = "red"))+ #NC
  theme_minimal()+
  ggtitle("HS model no nets (black), NC model no nets(red)")

#ivm_on <- 365*7
#year_period <- 3*365
#ivm_off <-  ivm_on + year_period
#


HS_NC_nets_EIR <- ggplot(df_3_epi, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_line(data = df_4_epi, aes(x = t, y = EIR_tot), col = "blue")+
  theme_minimal()+
  ggtitle("HS model (black) vs NC model (blue). Model with LLIN + IVM")



#abs and rel diff in prevalence: HS
#df_0_mean_prev <-  df_0_epi %>% # ivm only
#  filter(between(t, ivm_on, ivm_off)) %>%
#  summarise(mean_prev = mean(slide_prev0to5)) #0.6192437
#df_3_mean_prev <-  df_3_epi %>% # nets + ivm
#  filter(between(t, ivm_on, ivm_off)) %>%
#  summarise(mean_prev = mean(slide_prev0to5)) #0.3471302
#
#df_0_mean_prev- df_3_mean_prev # 0.2721134
#
#
#rel_prev_HS <- ((df_0_mean_prev - df_3_mean_prev)/df_0_mean_prev)*100 # 0.4394287
#
#df_2_mean_prev <-  df_2_epi %>% # ivm only
#  filter(between(t, ivm_on, ivm_off)) %>%
#  summarise(mean_prev = mean(slide_prev0to5)) #0.6172642
#
#df_4_mean_prev <-  df_4_epi %>% # ivm + nets
#  filter(between(t, ivm_on, ivm_off)) %>%
#  summarise(mean_prev = mean(slide_prev0to5)) #  0.3461478
#
#df_2_mean_prev - df_4_mean_prev #0.2711164
#
#rel_prev_NC = ((df_2_mean_prev - df_4_mean_prev)/df_2_mean_prev)*100 #43.92%
#
#HS_plots_prev <- ggplot(df_0_epi, aes(x = t, y = slide_prev0to5))+ #ivm only
#  geom_line()+
#  geom_line(data = df_3_epi, aes(x = t, y = slide_prev0to5), col = "red")+ # ivm + nets
#  geom_line(data = df_3a_epi, aes(x = t, y = slide_prev0to5), col = "blue")+ #nets only
#  ylim(0, 1)+
#  ggtitle("HS plots prev, black = ivm only, red = llin + ivm, blue = llin only")
#
#NC_plots_prev <- ggplot(df_2_epi, aes(x = t, y = slide_prev0to5))+
#  geom_line()+
#  geom_line(data = df_4_epi, aes(x = t, y = slide_prev0to5), col = "red")+ #ivm + nets
#  geom_line(data = df_2a_epi, aes(x = t, y = slide_prev0to5), col = "blue")+ # nets only
#  ylim(0, 1)+
#  ggtitle("NC plots prev, black = ivm only, red = llin + ivm, blue = llin only")
#
#plot_grid(HS_plots_prev, NC_plots_prev)
#
#HS_NC_nets_prev <- ggplot(df_3_epi, aes(x = t, y = slide_prev0to5))+
#  geom_line()+
#  geom_line(data = df_4_epi, aes(x = t, y = slide_prev0to5), col = "blue")+
#  theme_minimal()+
#  ggtitle("HS model (black) vs NC model (blue). Model with LLIN + IVM")+
#  ylim(0, 1)
#
#mod_comparison_plot <- plot_grid(HS_NC_nets_EIR, HS_NC_nets_prev)



#seems that the abs and rel EIRs do not change when we have a time-varying avhc and ivermectin vs fixed avhc and ivermectin
#this suggests that malariasimulation can work in a simple form and we do not need to account for the impact on delays in refeeding on excess mortality due to ivermectin
#mechanistically we can see (from the time-varying avhc) that LLINs change the proportion of mosquitoes killed from ivm (check this) but it doesn't seem to matter when it comes to examining the EIR or prevalence


#now check this across a range of LLIN coverages. NB that the mu_h will change with ivm coverage and bionomics, so keep these constant for now

#run HS and NC model for a range of LLIN coverages
#in NC model, we supply avhc_constant into avhc_in (inst/extdata/odin_endec_mu_h_constant_ivm_uptake.R)


HS_ITN_cov_loop <- function(itn_cov_in){
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = 100,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_parms1$ivm_cov_par, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    d_ITN0 = pyr_only_dn0,
    r_ITN0 = pyr_only_rn0
  )
  return(output)
}

itn_cov_vector <- seq(0, 0.8, 0.2)

HS_out_list <- lapply(itn_cov_vector, HS_ITN_cov_loop) #loop through all parameter values

res_HS_out_list <- lapply(HS_out_list, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
HS_out_df <- do.call(rbind,
                  sapply(1:length(itn_cov_vector), function(x){
                    as.data.frame(res_HS_out_list[[x]]) %>%
                      select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead) %>%
                      mutate(ref = x)
                  }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(HS_out_df, file = "data/llin_ivm_muh/HS_itn_cov_loop.csv", row.names = FALSE)


NC_ITN_cov_loop <- function(itn_cov_in){
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_endec_mu_h_constant_ivm_uptake.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = 100,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms3$haz, # hazard ratio for each off the eff_len number of days. SET TO 1
    ivm_cov_par = ivm_parms1$ivm_cov_par, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    avhc = avhc_constant,
    mu_h = 0.26,
    d_ITN0 = pyr_only_dn0,
    r_ITN0 = pyr_only_rn0
  )
  return(output)
}

NC_out_list <- lapply(itn_cov_vector, NC_ITN_cov_loop) #loop through all parameter values

res_NC_out_list <- lapply(NC_out_list, runfun) #put these values into the model

#now go through the res_out_list and save key parameters: av_mosq, LLIN cov, t, Sxtot, Extot, Ixtot
require(tidyverse)
NC_out_df <- do.call(rbind,
                     sapply(1:length(itn_cov_vector), function(x){
                       as.data.frame(res_NC_out_list[[x]]) %>%
                         select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, prop_killed_ivm, mv_dead, mvx_dead) %>%
                         mutate(ref = x)
                     }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(NC_out_df, file = "data/llin_ivm_muh/NC_itn_cov_loop.csv", row.names = FALSE)


  #does the change in mu exceed the change in avhc, which is why we don't see a difference between explicit model of refeeding delay vs not
#PLOT 2
#plot of avhc over time - this is the normalised value and what clearly affects the ivm uptake
#av_mosq and avhc under different LLIN coverages
out.df <- read.csv("data/out_df.csv", header = TRUE) #made in llin_ivm.R

avhc_plot <- ggplot(out.df, aes(x = t, y = avhc, col = as.factor(itn_cov)))+
  geom_line(linewidth = 1)+
  ylim(0, 1)


mu_plot <- ggplot(out.df, aes(x = t, y = mu, col = as.factor(itn_cov)))+
  geom_line(linewidth = 1)+
  ylim(0, 1)

#between 0 ITN cov and 0.8 ITN cov, what is absolute and rel diff in mu and avhc across full 10y

avhc_summary1 <- out.df %>%
  select(t, avhc, mu, itn_cov) %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc),
            mean_mu = mean(mu))

avhc_summary <- avhc_summary1 %>%
  mutate(abs_diff_avhc  = mean_avhc[1] - mean_avhc,
         rel_diff_avhc = (mean_avhc[1] - mean_avhc)/mean_avhc[1],
         abs_diff_mu = 0.132-mean_mu,
         rel_diff_mu = (0.132-mean_mu)/0.132)

abs_diff_params <- ggplot(avhc_summary, aes(x = itn_cov, y = abs_diff_avhc, col = as.factor(itn_cov)))+
  geom_point(shape = 2, size = 4)+
  geom_point(data = avhc_summary, aes(x = itn_cov, y = abs_diff_mu, fill = as.factor(itn_cov)), shape = 4, size = 4)+
  ylab("Absolute difference (triangle = avhc, cross = mu) \n compared to net cov 0")+
  theme_minimal()+
  ylim(-1.2, 1)

rel_diff_params <- ggplot(avhc_summary, aes(x = itn_cov, y = rel_diff_avhc, col = as.factor(itn_cov)))+
  geom_point(shape = 2, size = 4)+
  geom_point(data = avhc_summary, aes(x = itn_cov, y = rel_diff_mu, fill = as.factor(itn_cov)), shape = 4, size = 4)+
  ylab("Relative difference (triangle = avhc, cross = mu) \n compared to net cov 0")+
  theme_minimal()+
  ylim(-1.2, 1)

plot_grid(abs_diff_params, rel_diff_params)

ggplot(avhc_summary, aes(x = itn_cov, y = mean_avhc))+
  geom_line()+
  geom_line(data = avhc_summary, aes(x = itn_cov, y = mean_mu))+
  geom_line(col = "red")+
  ylim(0, 1)

out_df_UTN <- read.csv("data/llin_ivm_muh/out_df_UTN.csv", header = TRUE)

avhc_summary_UTN_1 <- out_df_UTN %>%
  select(t, avhc, mu, itn_cov) %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc),
            mean_mu = mean(mu))

avhc_summary_UTN <- avhc_summary_UTN_1 %>%
  mutate(abs_diff_avhc  = mean_avhc[1] - mean_avhc,
         rel_diff_avhc = (mean_avhc[1] - mean_avhc)/mean_avhc[1],
         abs_diff_mu = 0.132-mean_mu,
         rel_diff_mu = (0.132-mean_mu)/0.132)

abs_diff_params_UTN <- ggplot(avhc_summary_UTN, aes(x = itn_cov, y = abs_diff_avhc, col = as.factor(itn_cov)))+
  geom_point(shape = 2, size = 4)+
  geom_point(data = avhc_summary_UTN, aes(x = itn_cov, y = abs_diff_mu, fill = as.factor(itn_cov)), shape = 4, size = 4)+
  ylab("Absolute difference (triangle = avhc, cross = mu) \n with UTN compared to net cov 0")+
  theme_minimal()+
  ylim(-1.2, 1)

rel_diff_params_UTN <- ggplot(avhc_summary_UTN, aes(x = itn_cov, y = rel_diff_avhc, col = as.factor(itn_cov)))+
  geom_point(shape = 2, size = 4)+
  geom_point(data = avhc_summary_UTN, aes(x = itn_cov, y = rel_diff_mu, fill = as.factor(itn_cov)), shape = 4, size = 4)+
  ylab("Relative difference (triangle = avhc, cross = mu) \n with UTN compared to net cov 0")+
  theme_minimal()+
  ylim(-1.2, 1)

abs_rel_diff_net_eff <- plot_grid(abs_diff_params, abs_diff_params_UTN, rel_diff_params,  rel_diff_params_UTN,
          ncol = 2, nrow =2, labels = c("LLIN", "UTN"))

mu_avhc_plot <- ggplot(avhc_summary_UTN, aes(x = itn_cov, y = mean_avhc))+
  geom_line(linetype = "dashed",  col = "black")+
  geom_line(data = avhc_summary, aes(x = itn_cov, y = mean_avhc), col = "black")+
  geom_line(data = avhc_summary_UTN, aes(x = itn_cov, y = mean_mu), linetype = "dashed", col = "red")+
  geom_line(data = avhc_summary, aes(x = itn_cov, y = mean_mu), col = "red")+
  ggtitle("Mean value across 10y,dashed = UTN, solid = LLIN, red = mean_mu, black = mean_avhc")

#exploring different net types
path_nets <- "C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy"
IG2_nets <- file.path(path_nets, "pyrethroid_pyrrole_nets.csv")
pyr_only <- file.path(path_nets, "pyrethroid_only_nets.csv")
pbo_nets <- file.path(path_nets, "pyrethroid_pbo_nets.csv")

IG2_nets_df <- read.csv(IG2_nets, header = TRUE)
pyr_only_df <- read.csv(pyr_only, header = TRUE)
pbo_nets_df <- read.csv(pbo_nets, header = TRUE)


#modify the dfs so that the gamma for pbos and pyrrole are same as pyr only and only have the no resistance scenario
IG2_nets_df <- IG2_nets_df %>%
  mutate(gamman_lo10 = pyr_only_df$gamman_lo10,
         gamman_med = pyr_only_df$gamman_med,
         gamman_up90 = pyr_only_df$gamman_up90) %>%
  filter(resistance == 0)


pbo_nets_df <- pbo_nets_df %>%
  mutate(gamman_lo10 = pyr_only_df$gamman_lo10,
         gamman_med = pyr_only_df$gamman_med,
         gamman_up90 = pyr_only_df$gamman_up90) %>%
  filter(resistance == 0)


pyr_only_df <- pyr_only_df %>%
  filter(resistance == 0)


IG2_d_ITN0 <- IG2_nets_df$dn0_med

pyr_d_ITN0 <- pyr_only_df$dn0_med

pbo_d_ITN0 <- pbo_nets_df$dn0_med

#ivm_cov_vector <- seq(0, 1, 0.25)

#make an empty list
IG2_param_list <- list()
pyr_param_list <- list()
pbo_param_list <- list()

#use expand grid to get all combinations of the parameter values
#param_df <- expand.grid(d_ITN0 = d_ITN0_vector,
#                        r_ITN0 = r_ITN0_vector) #,
#ivm_cov = ivm_cov_vector)

itn_cov_vector <- seq(0, 0.8, 0.2)

#only use expand grid on the ivm_cov and dn0, because rn0 is correlated with dn0 so we add it in via left_join
IG2_param_df_crit <- expand.grid(dn0_med = IG2_d_ITN0, itn_cov = itn_cov_vector) #606 rows

IG2_param_df <- left_join(IG2_param_df_crit,
                          IG2_nets_df %>% dplyr::select(dn0_med, rn0_med),
                          by = c("dn0_med")) %>%
  mutate(net_type = "IG2") %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med)

#repeat for pyrethroid nets
pyr_param_df_crit <- expand.grid(dn0_med = pyr_d_ITN0, itn_cov = itn_cov_vector)
pyr_param_df <- left_join(pyr_param_df_crit,
                          pyr_only_df %>% dplyr::select(dn0_med, rn0_med),
                          by = c("dn0_med")) %>%
  mutate(net_type = "pyrethroid only") %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med)

#repeat for pbo nets
pbo_param_df_crit <- expand.grid(dn0_med = pbo_d_ITN0, itn_cov = itn_cov_vector)
pbo_param_df <- left_join(pbo_param_df_crit,
                          pbo_nets_df %>% dplyr::select(dn0_med, rn0_med),
                          by = c("dn0_med")) %>%
  mutate(net_type = "pyrethroid pbo") %>%
  rename(d_ITN0 = dn0_med, r_ITN0 = rn0_med)

#for now, no resistance, but could do filtration here of other net types etc
IG2_param_df_filt <- IG2_param_df

pyr_param_df_filt <- pyr_param_df

pbo_param_df_filt <- pbo_param_df


#conversions to numeric for each df
for(i in seq_len(nrow(IG2_param_df_filt))){
  IG2_param_list[[i]] <- as.numeric(IG2_param_df_filt[i,]) #convert to numeric so it's in the right form for runfun
}

for(i in seq_len(nrow(pyr_param_df_filt))){
  pyr_param_list[[i]] <- as.numeric(pyr_param_df_filt[i,]) #convert to numeric so it's in the right form for runfun
}

for(i in seq_len(nrow(pbo_param_df_filt))){
  pbo_param_list[[i]] <- as.numeric(pbo_param_df_filt[i,]) #convert to numeric so it's in the right form for runfun
}


#function for looping through net coverages and net types

#don't need to hardcode the half life as it is in the model_parameters file
#will need to put changes in here when add in resistance etc

create_itn_type_vec <- function(itn_type_ivm_param){
  d_ITN0_in <- itn_type_ivm_param[1]
  itn_cov_in <-itn_type_ivm_param[2]
  r_ITN0_in <- itn_type_ivm_param[3]
  #itn_half_life_in <- itn_type_ivm_param[4]*365
  output <- ivRmectin::create_r_model(
    odin_model_path = "inst/extdata/odin_model_endectocide.R",
    #num_int = 1,
    num_int = 2, # number of vector control (IRS and ITN) population groups
    ITN_IRS_on = itn_on,
    itn_cov = itn_cov_in,
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
    ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
    ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start ,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
  )
  return(output)
}

IG2_out_list <- lapply(IG2_param_list, create_itn_type_vec) #putting param list into the function to generate parameter set
pyr_out_list <- lapply(pyr_param_list, create_itn_type_vec)
pbo_out_list <- lapply(pbo_param_list, create_itn_type_vec)


#run it
IG2_res_out_list <- lapply(IG2_out_list, runfun)
pyr_res_out_list <- lapply(pyr_out_list, runfun)
pbo_res_out_list <- lapply(pbo_out_list, runfun)


IG2_out_df<- do.call(rbind,
                       sapply(1:(length(itn_cov_vector)), function(x){
                         as.data.frame(IG2_res_out_list[[x]]) %>%
                           select(t, mv, mu, avhc, EIR_tot, itn_cov, prop_killed_ivm) %>%
                           mutate(ref = x,
                                  net_type = "IG2")
                       }, simplify = F))

write.csv(IG2_out_df, file = "data/llin_ivm_muh/IG2_df_nores.csv", row.names = FALSE)

pyr_out_df <- do.call(rbind,
                        sapply(1:(length(itn_cov_vector)), function(x){
                          as.data.frame(pyr_res_out_list[[x]]) %>%
                            select(t, mv,mu, avhc, EIR_tot, itn_cov, prop_killed_ivm) %>%
                            mutate(ref = x,
                                   net_type = "pyr_only")
                        }, simplify = F))

write.csv(pyr_out_df, file = "data/llin_ivm_muh/pyr_df_nores.csv", row.names = FALSE)

pbo_out_df <- do.call(rbind,
                        sapply(1:(length(itn_cov_vector)), function(x){
                          as.data.frame(pbo_res_out_list[[x]]) %>%
                            select(t, mv,mu, avhc, EIR_tot, itn_cov, prop_killed_ivm) %>%
                            mutate(ref = x,
                                   net_type = "pbo")
                        }, simplify = F))

write.csv(pbo_out_df, file = "data/llin_ivm_muh/pbo_df_nores.csv", row.names = FALSE)

#look at avhc and mu in these net types

all_nets <- do.call("rbind",list(IG2_out_df, pbo_out_df, pyr_out_df))

all_nets_plot <- ggplot(all_nets, aes(x = t, y = EIR_tot, col = as.factor(itn_cov), linetype = net_type))+
  geom_line()

IG2_summary <- IG2_out_df %>%
  select(t, avhc, mu, itn_cov) %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc),
            mean_mu = mean(mu)) %>%
  mutate(abs_diff_avhc  = mean_avhc[1] - mean_avhc,
         rel_diff_avhc = (mean_avhc[1] - mean_avhc)/mean_avhc[1],
         abs_diff_mu = mean_mu[1]-mean_mu,
         rel_diff_mu = (mean_mu[1]-mean_mu)/mean_mu[1],
         net_type = "IG2")

pyr_summary <- pyr_out_df %>%
  select(t, avhc, mu, itn_cov) %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc),
            mean_mu = mean(mu)) %>%
  mutate(abs_diff_avhc  = mean_avhc[1] - mean_avhc,
         rel_diff_avhc = (mean_avhc[1] - mean_avhc)/mean_avhc[1],
         abs_diff_mu = mean_mu[1]-mean_mu,
         rel_diff_mu = (mean_mu[1]-mean_mu)/mean_mu[1],
         net_type = "pyrethroid")

pbo_summary <- pbo_out_df %>%
  select(t, avhc, mu, itn_cov) %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc),
            mean_mu = mean(mu)) %>%
  mutate(abs_diff_avhc  = mean_avhc[1] - mean_avhc,
         rel_diff_avhc = (mean_avhc[1] - mean_avhc)/mean_avhc[1],
         abs_diff_mu = mean_mu[1]-mean_mu,
         rel_diff_mu = (mean_mu[1]-mean_mu)/mean_mu[1],
         net_type = "pbo")

nets_summary <- do.call("rbind", list(IG2_summary, pyr_summary, pbo_summary))


abs_diff_net_summary <- ggplot(nets_summary, aes(x = itn_cov, y = abs_diff_avhc, linetype = as.factor(net_type)))+
  geom_line()+
  geom_line(data = nets_summary, aes(x = itn_cov, y = abs_diff_mu, linetype = as.factor(net_type)), col = "red")+
  geom_line()+
  ylim(-1, 1)+
  ylab("Abs diff in avhc (black) and mu (red) from net cov 0")


rel_diff_net_summary <- ggplot(nets_summary, aes(x = itn_cov, y = rel_diff_avhc, linetype = as.factor(net_type)))+
  geom_line()+
  geom_line(data = nets_summary, aes(x = itn_cov, y = rel_diff_mu, linetype = as.factor(net_type)), col = "red")+
  geom_line()+
  ylim(-1, 1)+
  ylab("Rel diff in avhc (black) and mu (red) from net cov 0")

plot_grid(abs_diff_net_summary, rel_diff_net_summary)




# looking at net distribution times: does this affect the reduction in the EIR from baseline and difference between models?

itn_on = 100

itn_distrib_times <- seq(itn_on, time_period, 365*3) #4 distributions in this period: d100 1195 2290 3385

#find the minimums

min_EIR_3a_1 <- df_3a_epi %>%
  filter(between(t, itn_distrib_times[1], itn_distrib_times[2]-1)) %>%
  slice_min(EIR_tot)

min_EIR_3a_2 <- df_3a_epi %>%
  filter(between(t, itn_distrib_times[2], itn_distrib_times[3]-1)) %>%
  slice_min(EIR_tot)

min_EIR_3a_3 <- df_3a_epi %>%
  filter(between(t, itn_distrib_times[3], itn_distrib_times[4]-1)) %>%
  slice_min(EIR_tot)
#rbind these
min_EIR_t_3a <- do.call("rbind", list(min_EIR_3a_1, min_EIR_3a_2, min_EIR_3a_3))
class(min_EIR_t_3a)


ggplot(df_3a_epi, aes( x = t, y = EIR_tot))+
  geom_point()+
  geom_point(data = itn_times_3a, aes( x = t, y = EIR_tot), col = "red")+
  geom_point(data = min_EIR_t_3a, aes(x = t, y = EIR_tot), col = "blue")


#what happens when you distribute IVM when nets have just been release (after the red points)

#scenario where you distribute ivm 10 days after an LLIN distribution

itn_distrib_round3 <- itn_distrib_times[3]

#start IVM in y6
ivm_parms2 <- ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = c(itn_distrib_round3+10, itn_distrib_round3+10+30, itn_distrib_round3+10+60),
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

#then run Hannah's models with LLINs (cov = 60%) and output the avhc
wh7 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms2$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  #eff_len = 60,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms2$ivm_cov_par,
                                  ivm_min_age = ivm_parms2$ivm_min_age,
                                  ivm_max_age = ivm_parms2$ivm_max_age,
                                  IVRM_start = ivm_parms2$IVRM_start)

#run Hannah's model
res7 <- runfun(wh7)
df_7 <- as.data.frame(res7)
write.csv(df_7, file = "data/llin_ivm_muh/df_7.csv", row.names = FALSE)

df_7_epi <- read.csv("data/llin_ivm_muh/df_7.csv", header = TRUE)

on <- (round((itn_distrib_round3+10)/365)-1)*365 #y5
off <- on + (3*365)

df7_mean_EIR <- df_7_epi %>% #ivm only
  filter(between(t, on, off)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 33.65784

ggplot(df_7_epi, aes(x = t, y = EIR_tot))+
  geom_point()

#same for the NC model

wh8 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_endec_mu_h_constant_ivm_uptake.R",
                                  #num_int = 1,
                                  num_int = 2,
                                  itn_cov = 0.6,
                                  ITN_IRS_on = 100,
                                  #het_brackets = 5,
                                  #age = init_age,
                                  init_EIR = 100,
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick",
                                  ttt = ivm_parms2$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  #eff_len = 60,
                                  haz = ivm_parms3$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  mu_h = 0.26,
                                  avhc = avhc_constant,
                                  IVRM_start = ivm_parms2$IVRM_start,
                                  d_ITN0 = pyr_only_dn0,
                                  r_ITN0 = pyr_only_rn0) #change the ivm on time

#run Hannah's model with mu_h of 0.26
res8 <- runfun(wh8)
df_8 <- as.data.frame(res8)
write.csv(df_8, file = "data/llin_ivm_muh/df_8.csv", row.names = FALSE)

df_8_epi <- read.csv("data/llin_ivm_muh/df_8.csv", header = TRUE)

df8_mean_EIR <- df_8_epi %>% #ivm only
  filter(between(t, on, off)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 33.62864


#even when you distribute IVM at the same time as a net distribution, the NC model can capture the EIRs given by Hannah's model
ggplot(df_7_epi, aes(x = t, y = EIR_tot))+
  geom_point()+
  geom_point(data = df_8_epi, aes(x = t, y = EIR_tot), col = "red")

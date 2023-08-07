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

#mv during the ivm distrib time
#then select mv during the ivermectin distribution time
#from d180 to d240+23
eff_len <- 23
ivm_on <- IVM_start[1]
ivm_off <- IVM_start[3]+eff_len
df_0_ivm <- df_0 %>%
  filter(between(t, ivm_on, ivm_off))

#then run the simplified version for different values of mu_h and see which values gives best fit to Hannah's mosquito densities
ivm_haz2 = rep(1, time_period) #not using hazards, so just set hazards to 1
ivm_parms3 <- ivm_fun(IVM_start_times = IVM_start,# time endectocide delivery occurs
                      #IVM_start = 180,
                      time_period = time_period,         # time period for the model to run over
                      hazard_profile = ivm_haz2[1:23], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                      #hazard_profile = ivm_haz2[1:730],
                      ivm_coverage = ivm_cov, # proportion of population receiving the endectocide
                      ivm_min_age = 5, # youngest age group receiving endectocide
                      ivm_max_age = 90) # oldest age group receiving endectocide

ivm_parms3$haz

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
                                  IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model with mu_h of 0.26
res2a <- runfun(wh2a)
df_2a <- as.data.frame(res2a)

#can see that mv broadly matches
mv_HS_NC <- ggplot(df_0, aes(x = t, y = mv))+ #HS
  geom_line()+
  geom_line(data = df_2, aes(x = t, y = mv, col = "red"), show.legend = FALSE)+ #NC
  theme_minimal()+
  ggtitle("mv w HS model (black) and NC model (red)")
ggsave(mv_HS_NC, file = "plots/llin_ivm_muh/mv_HS_NC.png")

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
                                  IVRM_start = ivm_parms1$IVRM_start)

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
                                  IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model
res3a <- runfun(wh3a)
df_3a <- as.data.frame(res3a)


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
                                  IVRM_start = ivm_parms1$IVRM_start)

#run Hannah's model with mu_h of 0.26
res4 <- runfun(wh4)
df_4 <- as.data.frame(res4)



#df_0 : Hannah's model, no nets, ivm
df_0_epi <- df_0 %>%
  select(t, EIR_tot, slide_prev0to5)
write.csv(df_0_epi, file = "data/llin_ivm_muh/df_0_epi.csv", row.names = FALSE)

#df_3: Hannah's model, nets, ivm
df_3_epi <- df_3 %>%
  select(t, EIR_tot, slide_prev0to5)
write.csv(df_3_epi, file = "data/llin_ivm_muh/df_3_epi.csv", row.names = FALSE)

#df_3a: Hannah's model, nets, no ivm
df_3a_epi <- df_3a %>%
  select(t, EIR_tot, slide_prev0to5)
write.csv(df_3a_epi, file = "data/llin_ivm_muh/df_3a_epi.csv", row.names = FALSE)

#df2: Nilani's model, no nets, ivm
df_2_epi <- df_2 %>%
  select(t, EIR_tot, slide_prev0to5)
write.csv(df_2_epi, file = "data/llin_ivm_muh/df_2_epi.csv", row.names = FALSE)

#df2a: Nilani's model, nets only
df_2a_epi <- df_2a %>%
  select(t, EIR_tot, slide_prev0to5)
write.csv(df_2a_epi, file = "data/llin_ivm_muh/df_2a_epi.csv", row.names = FALSE)


#df4: Nilani's model, nets, ivm
df_4_epi <- df_4 %>%
  select(t, EIR_tot, slide_prev0to5)
write.csv(df_4_epi, file = "data/llin_ivm_muh/df_4_epi.csv", row.names = FALSE)

#READ IN THE EPI FILES HERE

#HS models

#HS, IVM only
df_0_epi <- read.csv("data/llin_ivm_muh/df_0_epi.csv", header = TRUE)

#HS, IVM + NETS
df_3_epi <- read.csv("data/llin_ivm_muh/df_3_epi.csv", header = TRUE)


#HS, NETS only
df_3a_epi <- read.csv("data/llin_ivm_muh/df_3a_epi.csv", header = TRUE)


#NC models

#NC, IVM only
df_2_epi <- read.csv("data/llin_ivm_muh/df_2_epi.csv", header = TRUE)


#NC, IVM + NETS
df_4_epi <- read.csv("data/llin_ivm_muh/df_4_epi.csv", header = TRUE)


#NC, NETS only
df_2a_epi <- read.csv("data/llin_ivm_muh/df_2a_epi.csv", header = TRUE)

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

year_ivm_before_on <- 365*7
year_period <- 3*365
y3_end_ivm <-  year_ivm_on + year_period

#absolute difference in EIR between df_0 and df_3
df_0_mean_EIR = df_0_epi %>% #ivm only
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 69.41213


df_3_mean_EIR = df_3_epi %>% #nets + ivm
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_EIR = mean(EIR_tot))  # 31.58684. nets + ivm

df_0_mean_EIR - df_3_mean_EIR #37.82529


df_2_mean_EIR = df_2_epi %>% #ivm only
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 69.01957. ivm only

df_4_mean_EIR = df_4_epi %>% # nets + ivm
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_EIR = mean(EIR_tot)) # 31.40277. nets + ivm

df_2_mean_EIR - df_4_mean_EIR # 37.6168

#rel diff HS
rel_EIR_HS <- ((df_0_mean_EIR-df_3_mean_EIR)/df_0_mean_EIR)*100 #54.49378

#rel diff NC
rel_EIR_NC <- ((df_2_mean_EIR-df_4_mean_EIR)/df_2_mean_EIR)*100 # 54.50164


HS_plots_EIR <- ggplot(df_0_epi, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_line(data = df_3_epi, aes(x = t, y = EIR_tot), col = "red")+
  geom_line(data = df_3a_epi, aes(x = t, y = EIR_tot), col = "blue")+
  ggtitle("EIR w HS model, ivm only (black), LLIN + IVM (red) and nets only (blue)")

NC_plots_EIR <- ggplot(df_2_epi, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_line(data = df_4_epi, aes(x = t, y = EIR_tot), col = "red")+
  geom_line(data = df_2a_epi, aes(x = t, y = EIR_tot), col = "blue")+
  ggtitle("EIR w NC model, iivm only (black), LLIN + IVM (red) and nets only (blue)")

require(cowplot)
plot_grid(HS_plots_EIR, NC_plots_EIR)

HS_NC_nets_EIR <- ggplot(df_3_epi, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_line(data = df_4_epi, aes(x = t, y = EIR_tot), col = "blue")+
  theme_minimal()+
  ggtitle("HS model (black) vs NC model (blue). Model with LLIN + IVM")



#abs and rel diff in prevalence: HS
df_0_mean_prev <-  df_0_epi %>% # ivm only
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_prev = mean(slide_prev0to5)) #0.6192437
df_3_mean_prev <-  df_3_epi %>% # nets + ivm
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_prev = mean(slide_prev0to5)) #0.3471302

df_0_mean_prev- df_3_mean_prev # 0.2721134


rel_prev_HS <- ((df_0_mean_prev - df_3_mean_prev)/df_0_mean_prev)*100 # 0.4394287

df_2_mean_prev <-  df_2_epi %>% # ivm only
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_prev = mean(slide_prev0to5)) #0.6172642

df_4_mean_prev <-  df_4_epi %>% # ivm + nets
  filter(between(t, year_ivm_before_on, y3_end_ivm)) %>%
  summarise(mean_prev = mean(slide_prev0to5)) #  0.3461478

df_2_mean_prev - df_4_mean_prev #0.2711164

rel_prev_NC = ((df_2_mean_prev - df_4_mean_prev)/df_2_mean_prev)*100 #43.92%

HS_plots_prev <- ggplot(df_0_epi, aes(x = t, y = slide_prev0to5))+ #ivm only
  geom_line()+
  geom_line(data = df_3_epi, aes(x = t, y = slide_prev0to5), col = "red")+ # ivm + nets
  geom_line(data = df_3a_epi, aes(x = t, y = slide_prev0to5), col = "blue")+ #nets only
  ylim(0, 1)+
  ggtitle("HS plots prev, black = ivm only, red = llin + ivm, blue = llin only")

NC_plots_prev <- ggplot(df_2_epi, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  geom_line(data = df_4_epi, aes(x = t, y = slide_prev0to5), col = "red")+ #ivm + nets
  geom_line(data = df_2a_epi, aes(x = t, y = slide_prev0to5), col = "blue")+ # nets only
  ylim(0, 1)+
  ggtitle("NC plots prev, black = ivm only, red = llin + ivm, blue = llin only")

plot_grid(HS_plots_prev, NC_plots_prev)

HS_NC_nets_prev <- ggplot(df_3_epi, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  geom_line(data = df_4_epi, aes(x = t, y = slide_prev0to5), col = "blue")+
  theme_minimal()+
  ggtitle("HS model (black) vs NC model (blue). Model with LLIN + IVM")+
  ylim(0, 1)

mod_comparison_plot <- plot_grid(HS_NC_nets_EIR, HS_NC_nets_prev)



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
    IVRM_start = ivm_parms1$IVRM_start
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
                      select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5) %>%
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
    mu_h = 0.26
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
                         select(t, mv, avhc, itn_cov, EIR_tot, slide_prev0to5) %>%
                         mutate(ref = x)
                     }, simplify = F))

#go from wide to long
#out_df_long <- gather(out_df, parameter, value, itn_cov:ref, factor_key = TRUE)

write.csv(NC_out_df, file = "data/llin_ivm_muh/NC_itn_cov_loop.csv", row.names = FALSE)

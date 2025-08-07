#here, I am distributing the endectocide against before the killing effect has waned
#the endectocide has an additive killing-effect of 0.01 and kills mosquitoes for 10 days
#i am going to distribute it to the population on d11, d14 and d17

#compare this to something that is:

#i) just distributed once
#ii) distributed with bigger gaps so that the efficacy has waned by the time of the next distribution: d11, d24, d37

#script to compare different model types
devtools::load_all()
require(tidyverse)
#how much does mu_h change across different values of Q0 and ivm_cov?

# Provide a value of the annual EIR for this model run
#init_EIR_vec <- c(2, 25, 100) #low - 2, moderate - 15, high - 120 --> Ellie: low = 2, med = 25, high = 100

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years, turn ivermectin on when nets are 6m, 1y, 2.5yo
time_period <- 100

#ivm on when nets are 6 months old
IVM_start <- c(11, 40, 90)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

#itn_on <- 100 #introduce nets 100 days into simulation


runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#here the hazard profile is 1 everyday

#set up the different start times (coverage is redundant here, we are just operating a switch)

ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution once every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

mod1 <- ivRmectin:::create_r_model(
  odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
  num_int = 1,
  #num_int = 2,
  #ITN_IRS_on = 100,
  #itn_cov = 0.75,
  #het_brackets = 5,
  #age = init_age,
  init_EIR = 25,
  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
  #admin2 = "Fatick",
  ttt = ivm_parms1$ttt,
  eff_len = ivm_parms1$eff_len,
  haz = ivm_parms1$haz,
  ivm_cov_par = ivm_parms1$ivm_cov_par,
  ivm_min_age = ivm_parms1$ivm_min_age,
  ivm_max_age = ivm_parms1$ivm_max_age,
  IVRM_start = ivm_parms1$IVRM_start,
  Q0 = 0.9)

mod1_df <- as.data.frame(runfun(mod1)) %>%
  select(t, mu, mv, mvx_dead, EIR_tot, slide_prev0to5, IVRM_sr)

ggplot(mod1_df, aes(x = t, y = mv))+
  geom_line()

mod1_endec_killed <- mod1_df %>%
  filter(between(t, IVM_start[1], IVM_start[1]+10)) %>%
  select(mv_dead) %>%
  summarise(tot_mv_dead = sum(mv_dead)) #have 90,786.99 mosquitoes killed in this period

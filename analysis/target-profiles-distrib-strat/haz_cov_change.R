require(tidyverse)
ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")

ivm_haz <- ivm_haz[1:23,]

ggplot(ivm_haz, aes(x = Day, y = IVM_300_3_HS))+
  geom_point()+
  geom_hline(yintercept = 1, col = "red", linetype = "dashed")

#two distributions, 80% target coverage gets split to two evenly sized distributions
#compute the eligible population
ivm_cov_par <- 0.8
ivm_min_age <- 5
ivm_max_age <- 90
ivm_cov <- ivm_cov_par*(exp(-ivm_min_age/21) - exp(-ivm_max_age/21))
#so ivm_cov = 0.619..


#say they treat the first population on day 1,2,3, then move to the next population on day 5
ivm_cov_haz <- ivm_haz %>%
  select(Day, IVM_300_3_HS) %>%
  rename(haz_300 = IVM_300_3_HS) %>%
  mutate(haz_300_A = haz_300*ivm_cov) #is this correct
#write.csv(ivm_haz, file = "IVM_derivation/ivm_hazards_change.csv")

#then compare runs with 80% coverage and this
runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

time_period <- 365
ivm_parms1 <- ivm_fun(IVM_start_times = 100,            # time endectocide delivery occurs
                     time_period = time_period,         # time period for the model to run over
                     hazard_profile = ivm_cov_haz$haz_300, # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                     ivm_coverage = 0.8, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide

ivm_parms2 <- ivm_fun(IVM_start_times = 100,            # time endectocide delivery occurs
                      time_period = time_period,         # time period for the model to run over
                      hazard_profile = ivm_cov_haz$haz_300_A, # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                      ivm_coverage = 1, # proportion of population receiving the endectocide
                      ivm_min_age = 5, # youngest age group receiving endectocide
                      ivm_max_age = 90) # oldest age group receiving endectocide

mod1 <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                num_int = 1, # number of vector control (IRS and ITN) population groups
                                #het_brackets = 5, # number of heterogeneous biting categories
                                #age = init_age, # the different age classes to be ran within the model
                                init_EIR = 100, # the Entomological Innoculation Rate
                                #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
                                eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
                                haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
                                ivm_cov_par = ivm_parms1$ivm_cov_par, # proportion of popuulation receiving the endectocide
                                ivm_min_age = ivm_parms1$ivm_min_age, # youngest age group receiving endectocide
                                ivm_max_age = ivm_parms1$ivm_max_age, # oldest age group receiving endectocide
                                IVRM_start = ivm_parms1$IVRM_start) # model specific parameter to control timing of endectocide delivery
mod1_out <- runfun(mod1)

mod2 <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 1, # number of vector control (IRS and ITN) population groups
                                  #het_brackets = 5, # number of heterogeneous biting categories
                                  #age = init_age, # the different age classes to be ran within the model
                                  init_EIR = 100, # the Entomological Innoculation Rate
                                  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  ttt = ivm_parms2$ttt, # model specific parameter to control timing of endectocide delivery
                                  eff_len = ivm_parms2$eff_len, # number of days after receiving endectocide that HR is higher
                                  haz = ivm_parms2$haz, # hazard ratio for each off the eff_len number of days
                                  ivm_cov_par = ivm_parms2$ivm_cov_par, # proportion of popuulation receiving the endectocide
                                  ivm_min_age = ivm_parms2$ivm_min_age, # youngest age group receiving endectocide
                                  ivm_max_age = ivm_parms2$ivm_max_age, # oldest age group receiving endectocide
                                  IVRM_start = ivm_parms2$IVRM_start) # model specific parameter to control timing of endectocide delivery
mod1_out <- runfun(mod1)
mod2_out <- runfun(mod2)

mod1_out <- as.data.frame(mod1_out) %>%
  mutate(model = "hazards_original")

mod2_out <- as.data.frame(mod2_out) %>%
  mutate(model = "hazards_scaled")

mods <- rbind(mod1_out, mod2_out)

ggplot(mods, aes(x = t, y = EIR_tot, col = as.factor(model)))+
  geom_line()

ggplot(mods, aes(x = t, y = mv, col = as.factor(model)))+
  geom_line()

ggplot(mods, aes(x = t, y = slide_prev0to5, col = as.factor(model)))+
  geom_line()+
  ylim(0,1)

#not sure how to deal with the coverage parameter



View(admin_units_seasonal)
setwd("~/OneDrive - LSTM/ivRmectin")
# Loading the ivRmectin package
devtools::load_all()


# Create a vector of age categories for the model
init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 3

#low - 2, moderate - 15, high - 120


# Provide the length of time (in days) that you want to run the model for
time_period <- 365 * 10# run model for 10 years

# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

# Load ivermectin hazardd
sim_ave_HR <- read.csv("/Users/anna.trett/OneDrive - LSTM/CHIC 599/NTBC_Vs_IVM_Summer2020/PKProfiles/01-Data/sim_ave_HR.csv", header=TRUE)

# Running the ivm_fun function to generate the extra endectocide specific parameters that
# you have to pass to the model
ivm_parms = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                     time_period = time_period, # time period for the model to run over
                     hazard_profile = sim_ave_HR$ave_HR_1[1:28], # dummy hazard profile - must be vector
                     ivm_coverage = 0.8, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide

# Creates the odin model with all the required parameters - it is then ready to run
# although this isn't the part where the model is actually run - that's below.
wh <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                num_int = 4, # number of vector control (IRS and ITN) population groups
                                het_brackets = 5, # number of heterogeneous biting categories
                                age = init_age, # the different age classes to be ran within the model
                                init_EIR = init_EIR, # the Entomological Innoculation Rate
                                country = NULL, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                admin2 = NULL, # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                ttt = ivm_parms$ttt, # model specific parameter to control timing of endectocide delivery
                                eff_len = ivm_parms$eff_len, # number of days after receiving endectocide that HR is higher
                                haz = ivm_parms$haz, # hazard ratio for each off the eff_len number of days
                                ivm_cov_par = ivm_parms$ivm_cov_par, # proportion of popuulation receiving the endectocide
                                ivm_min_age = ivm_parms$ivm_min_age, # youngest age group receiving endectocide
                                ivm_max_age = ivm_parms$ivm_max_age, # oldest age group receiving endectocide
                                IVRM_start = ivm_parms$IVRM_start) # model specific parameter to control timing of endectocide delivery

# This function runs the odin model with all the parameters that you previously generated with
# the ivRmectimn::create_r_model call
runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

runfun(wh)

sim_ave_HR -> sim_ave_HR %>%
    mutate(sim_ave_HR$HR=4)
View(sim_ave_HR)
sim_ave_HR$HR=4
# Running the ivm_fun function to generate the extra endectocide specific parameters that
# you have to pass to the model
ivm_parms = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                     time_period = time_period, # time period for the model to run over
                     hazard_profile = sim_ave_HR$ave_HR_1[1:23], # dummy hazard profile - must be vector
                     ivm_coverage = 0.8, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide

# Creates the odin model with all the required parameters - it is then ready to run
# although this isn't the part where the model is actually run - that's below.
wh <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                num_int = 4, # number of vector control (IRS and ITN) population groups
                                het_brackets = 5, # number of heterogeneous biting categories
                                age = init_age, # the different age classes to be ran within the model
                                init_EIR = init_EIR, # the Entomological Innoculation Rate
                                country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                ttt = ivm_parms$ttt, # model specific parameter to control timing of endectocide delivery
                                eff_len = ivm_parms$eff_len, # number of days after receiving endectocide that HR is higher
                                haz = ivm_parms$haz, # hazard ratio for each off the eff_len number of days
                                ivm_cov_par = ivm_parms$ivm_cov_par, # proportion of popuulation receiving the endectocide
                                ivm_min_age = ivm_parms$ivm_min_age, # youngest age group receiving endectocide
                                ivm_max_age = ivm_parms$ivm_max_age, # oldest age group receiving endectocide
                                IVRM_start = ivm_parms$IVRM_start) # model specific parameter to control timing of endectocide delivery


#graphical output

#COmparison of age coverage between IVM maximum dose
#Adult incidence and prevalence outcomes

#highly seasonal
#low burden
Incidence = Params_Func (eir=3, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_1[1:23], cov2=0.8, min2=0, max2=90,
                         ylim1=c(0, 1500), ylim2= c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))

#moderate bburden
Incidence = Params_Func (eir=14,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_1[1:23], cov2=0.8, min2=0, max2=90,
                         ylim1=c(0, 4000), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))

#high burden
Incidence = Params_Func (eir=500,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_1[1:23], cov2=0.8, min2=0,  max2=90,
                         ylim1=c(0, 9000), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))



#Perrenial seasonal
#low burden
Incidence = Params_Func (eir=3,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_1[1:23], cov2=0.8, min2=0, max2=90,
                         ylim1=c(0, 1500), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))
#moderate burden
Incidence = Params_Func (eir=14,   start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_1[1:23], cov2=0.8, min2=0, max2=90,
                         ylim1=c(0, 1500), ylim2=c(0,70), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))
#high burden
Incidence = Params_Func (eir=500,   start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_1[1:23], cov2=0.8, min2=0, max2=90,
                         ylim1=c(0, 2500), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "IVM (0-5yrs)"))
################################################
#Comparing NTBC and IVM maximum doses, same age category coverage (5-90yrs)

#highly seasonal
#low burden
Incidence = Params_Func (eir=3,   start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=5, max2=90,
                         ylim1=c(0, 1500), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))


#moderate bburden
Incidence = Params_Func (eir=14,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=5, max2=90,
                         ylim1=c(0, 1500), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", " NTBC (5-90yrs)"))

#high burden
Incidence = Params_Func (eir=500,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=5, max2=90,
                         ylim=c(0, 8000), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))


#Perrenial seasonal
#low burden
Incidence = Params_Func (eir=3,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=5, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))

#moderate burden
Incidence = Params_Func (eir=14, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=5, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))

#high burden
Incidence = Params_Func (eir=500, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=5, max2=90,
                         ylim=c(0, 2500), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (5-90yrs)"))

################################################
#Comparing NTBC and IVM maximum doses, different age category coverage (5-90yrs)

#highly seasonal
#low burden
Incidence = Params_Func (eir=3, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-90yrs)"))

#moderate bburden
Incidence = Params_Func (eir=14,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-90yrs)"))

#high burden
Incidence = Params_Func (eir=500,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 9000), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-90yrs)"))


#Perrenial seasonal
#low burden
Incidence = Params_Func (eir=3,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1000), ylim2=c(0,80), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-90yrs)"))

#moderate burden
Incidence = Params_Func (eir=14, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-90yrs)"))

#high burden
Incidence = Params_Func (eir=500, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=5, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 3000), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "NTBC (0-90yrs)"))


################################################
#Comparing NTBC, different dosing regimen for NTBC

#highly seasonal
#low burden
Incidence = Params_Func (eir=3, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No endectocide", "NTBC (1x70mg)", "NTBC (3x70mg)"))

#moderate bburden
Incidence = Params_Func (eir=14,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80),labels=c("No endectocide", "NTBC (1x70mg)", "NTBC (3x70mg)"))

#high burden
Incidence = Params_Func (eir=500,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 8500), ylim2=c(0,100),labels=c("No endectocide", "NTBC (1x70mg)", "NTBC (3x70mg)"))


#Perrenial seasonal
#low burden
Incidence = Params_Func (eir=3,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1000), ylim2=c(0,70), labels=c("No endectocide", "NTBC (1x70mg)", "NTBC (3x70mg)"))

#moderate burden
Incidence = Params_Func (eir=14, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No endectocide", "NTBC (1x70mg)", "NTBC (3x70mg)"))

#high burden
Incidence = Params_Func (eir=500, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 3500), ylim2=c(0,100), labels=c("No endectocide", "NTBC (1x70mg)", "NTBC (3x70mg)"))




################################################
#Comparing NTBC, different dosing regimen for NTBC - Comparing timing of drug administration

#highly seasonal
#low burden
Incidence = Params_Func (eir=3, start1=c(3650, 3680, 4100), start2=c(3650, 3680, 4100), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1000), ylim2=c(0,80), labels=c("Without endectocide", "NTBC (1x70mg) (0-90yrs)", "NTBC (3x70mg) (0-90yrs)"))

#moderate bburden
Incidence = Params_Func (eir=14,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), labels=c("No endectocide", "NTBC (1x70mg) (0-90yrs)", "NTBC (3x70mg) (0-90yrs)"))

#high burden
Incidence = Params_Func (eir=500,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), labels=c("No endectocide", "NTBC (1x70mg) (0-90yrs)", "NTBC (3x70mg) (0-90yrs)"))

#Perrenial seasonal
#low burden
Incidence = Params_Func (eir=3,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), labels=c("No endectocide", "NTBC (1x70mg) (0-90yrs)", "NTBC (3x70mg) (0-90yrs)"))

#moderate burden
Incidence = Params_Func (eir=14, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), labels=c("No endectocide", "NTBC (1x70mg) (0-90yrs)", "NTBC (3x70mg) (0-90yrs)"))

#high burden
Incidence = Params_Func (eir=500, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_3[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_3[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$ave_HR_4[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 3000), labels=c("No endectocide", "NTBC (1x70mg) (0-90yrs)", "NTBC (3x70mg) (0-90yrs)"))



################
#TCP-6 profile versus IVM

#highly seasonal
#low burden
Incidence = Params_Func (eir=3, start1=c(3650, 3680, 4100), start2=c(3650, 3680, 4100), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No  endectocide", "IVM (5-90yrs)", "TCP-6 (0-90yrs)"))

#moderate bburden
Incidence = Params_Func (eir=14,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 4000), ylim2=c(0,80), labels=c("No endectocide", "IVM (5-90yrs)", "TCP-6 (0-90yrs)"))

#high burden
Incidence = Params_Func (eir=500,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 9000), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "TCP-6 (0-90yrs)"))


#Perennial seasonal
#low burden
Incidence = Params_Func (eir=3,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1000), ylim2=c(0,80), labels=c("No endectocide", "IVM (5-90yrs)", "TCP-6 (0-90yrs)"))

#moderate burden
Incidence = Params_Func (eir=14, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 2500), ylim2=c(0,80), labels=c("No endectocide", "IVM (5-90yrs)", "TCP-6 (0-90yrs)"))

#high burden
Incidence = Params_Func (eir=500, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_1[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_1[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 3500), ylim2=c(0,100), labels=c("No endectocide", "IVM (5-90yrs)", "TCP-6 (0-90yrs)"))


############################
#TCP-6 versus NTBC high dose
#highly seasonal
#low burden
Incidence = Params_Func (eir=3, start1=c(3650, 3680, 4100), start2=c(3650, 3680, 4100), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_4[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_4[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No endectocide", "NTBC (3x70mg)", "TCP-6"))

#moderate bburden
Incidence = Params_Func (eir=14,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_4[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_4[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 4000), ylim2=c(0,80), labels=c("No endectocide", "NTBC (3x70mg)", "TCP-6"))

#high burden
Incidence = Params_Func (eir=500,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Senegal", Region="Fatick",
                         HR_sim0=sim_ave_HR$ave_HR_4[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_4[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 8000), ylim2=c(0,100), labels=c("No endectocide", "NTBC (3x70mg)", "TCP-6 "))

#Perrenial seasonal
#low burden
Incidence = Params_Func (eir=3,  start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_4[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_4[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1000), ylim2=c(0,80), labels=c("No endectocide", "NTBC (3x70mg)", "TCP-6"))

#moderate burden
Incidence = Params_Func (eir=14, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_4[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_4[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 1500), ylim2=c(0,80), labels=c("No endectocide", "NTBC (1x70mg) (0-90yrs)", "TCP-6 (0-90yrs)"))

#high burden
Incidence = Params_Func (eir=500, start1=c(2250, 2280, 2310), start2=c(2250, 2280, 2310), Country="Democratic Republic of the Congo", Region="Equateur",
                         HR_sim0=sim_ave_HR$ave_HR_4[1:23],cov0=0.8, min0=5, max0=90,
                         HR_sim1=sim_ave_HR$ave_HR_4[1:23], cov1=0.8, min1=0, max1=90,
                         HR_sim2=sim_ave_HR$HR[1:23], cov2=0.8, min2=0, max2=90,
                         ylim=c(0, 2500), ylim2=c(0,100), labels=c("No endectocide", "NTBC (1x70mg) (0-90yrs)", "TCP-6 (0-90yrs)"))











###############################################
Params_Func = function (eir, start1, start2, Country, Region,
                        HR_sim0, cov0, min0, max0,
                        HR_sim1, cov1, min1, max1,
                        HR_sim2, cov2, min2, max2,
                        ylim1, ylim2, labels) {
  #, HR_sim2, cov2, min2, max2) {

# no ivermectin being used
  ivm_parms0 = ivm_fun(IVM_start_times = 10000,  # time endectocide delivery occurs
                     time_period = time_period, # time period for the model to run over
                     hazard_profile = HR_sim0, # dummy hazard profile - must be vector
                     ivm_coverage = cov0, # proportion of population receiving the endectocide
                     ivm_min_age = min0, # youngest age group receiving endectocide
                     ivm_max_age = max0) # oldest age group receiving endectocide


wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = eir,
                                  country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = Region,
                                  ttt = ivm_parms0$ttt,
                                  eff_len = ivm_parms0$eff_len,
                                  haz = ivm_parms0$haz,
                                  ivm_cov_par = ivm_parms0$ivm_cov_par,
                                  ivm_min_age = ivm_parms0$ivm_min_age,
                                  ivm_max_age = ivm_parms0$ivm_max_age,
                                  IVRM_start = ivm_parms0$IVRM_start)

res <- runfun(wh)


# ivermectin being delivered at 3 time points, with HR of 2, uniformly lasting for 10 days
ivm_parms1 = ivm_fun(IVM_start_times = start1,
                     time_period = time_period,
                     hazard_profile = HR_sim1,
                     ivm_coverage=cov1,
                     ivm_min_age=min1,
                     ivm_max_age = max1)


wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = eir,
                                  country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = Region,
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start)


ivm_parms2 = ivm_fun(IVM_start_times = start2,
                     time_period = time_period,
                     hazard_profile = HR_sim2,
                     ivm_coverage=cov2,
                     ivm_min_age=min2,
                     ivm_max_age = max2)

wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = eir,
                                  country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = Region,
                                  ttt = ivm_parms2$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms2$ivm_cov_par,
                                  ivm_min_age = ivm_parms2$ivm_min_age,
                                  ivm_max_age = ivm_parms2$ivm_max_age,
                                  IVRM_start = ivm_parms2$IVRM_start)

# Running each of the models wh0 (no ivermectin), wh1 (ivermectin, HR = 2 for 10 days) and
# wh2 (ivermectin, HR = 2 for 28 days)
res0 <- runfun(wh0)
res1 <- runfun(wh1)
res2 <- runfun(wh2)

# Plotting the results
cols <- c("grey40", "deeppink2", "deepskyblue3")
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))

# Clinical Incidence
plot = plot(res0$t/365, res0$clin_inc0to5*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
          lwd = 2, col = cols[1], xlim = c(5.7, 8), ylim = ylim1, las = 1)
lines(res1$t/365, res1$clin_inc0to5*1000*365, lwd = 2, col = cols[2])
lines(res2$t/365, res2$clin_inc0to5*1000*365, lwd = 2, col = cols[3])
arrows(start1/365, -50, start2/365, 10, length = 0.05, lwd = 1, col = "goldenrod2")

legend("topright", labels,
       col = cols, lwd=2, bty="n", cex=0.8)

plot(res0$t/365, res0$slide_prev0to5*100, type="l", ylab="Slide prevalence (%)", xlab="Year", lwd=3, col=cols[1],
     xlim = c(5.7, 8), ylim = ylim2, las=1)
lines(res1$t/365, res1$slide_prev0to5*100, lwd=2, col=cols[2])
lines(res2$t/365, res2$slide_prev0to5*100, lwd=2, col=cols[3])
arrows(start1/365, -50, start2/365, 1, length = 0.05, lwd = 1, col = "goldenrod2")

legend("topright", labels,
       col = cols, lwd=2, bty="n", cex=0.8)



return (plot)

}













ivm_parms2 = ivm_fun(IVM_start_times = c(2250, 2280, 2310),
                     time_period = time_period,
                     hazard_profile = HR_sim2,
                     ivm_coverage=cov2,
                     ivm_min_age=min2,
                     ivm_max_age = max2)

wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = eir,
                                  country = Country, # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = Region,
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms2$ivm_cov_par,
                                  ivm_min_age = ivm_parms2$ivm_min_age,
                                  ivm_max_age = ivm_parms2$ivm_max_age,
                                  IVRM_start = ivm_parms2$IVRM_start)

# Running each of the models wh0 (no ivermectin), wh1 (ivermectin, HR = 2 for 10 days) and
# wh2 (ivermectin, HR = 2 for 28 days)
res0 <- runfun(wh0)
res1 <- runfun(wh1)
res2 <- runfun(wh2)

# Plotting the results
cols <- c("grey40", "deeppink2", "deepskyblue3")
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))

# Clinical Incidence
plot=plot(res0$t/365, res0$clin_inc0to80*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 2, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$clin_inc0to80*1000*365, lwd = 2, col = cols[2])
lines(res2$t/365, res2$clin_inc0to80*1000*365, lwd = 2, col = cols[3])
arrows(c(2250, 2280, 2310)/365, -50, c(2250, 2280, 2310)/365, 10, length = 0.15, lwd = 1, col = "goldenrod2")


plot(res0$t/365, res0$slide_prev0to80*100, type="l", ylab="Slide prevalence (%)", xlab="Year", lwd=3, col=cols[1],
     xlim = c(5.7, 8), ylim = c(0, 70), las=1)
lines(res1$t/365, res1$slide_prev0to80*100, lwd=3, col=cols[2])
lines(res2$t/365, res2$slide_prev0to80*100, lwd=3, col=cols[3])
arrows(c(2250, 2280, 2310)/365, -50, c(2250, 2280, 2310)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")

legend("topright", c("No endectocide", "IVM endectocide", "NTBC endectocide"),
       col = cols, lwd=3, bty="n", cex=0.8)
return (plot)

}

#########################################

ivm_parms0 = ivm_fun(IVM_start_times = 10000,
                     time_period = time_period,
                     hazard_profile = sim_ave_HR$ave_HR_1[1:28],
                     ivm_coverage=0.8,
                     ivm_min_age=0,
                     ivm_max_age = 90)

wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms0$ivm_cov_par,
                                  ivm_min_age = ivm_parms0$ivm_min_age,
                                  ivm_max_age = ivm_parms0$ivm_max_age,
                                  IVRM_start = ivm_parms0$IVRM_start)




ivm_parms1 = ivm_fun(IVM_start_times = c(3650, 3680, 4100),
                     time_period = time_period,
                     hazard_profile = sim_ave_HR$ave_HR_1[1:28],
                     ivm_coverage=0.8,
                     ivm_min_age=0,
                     ivm_max_age = 90)
ivm_parms1
wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start)



ivm_parms2 = ivm_fun(IVM_start_times = c(2250, 2280, 2310),
                     time_period = time_period,
                     hazard_profile = sim_ave_HR$ave_HR_1[1:28],
                     ivm_coverage=0.8,
                     ivm_min_age=0,
                     ivm_max_age = 90)

wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms2$eff_len,
                                  haz = ivm_parms2$haz,
                                  ivm_cov_par = ivm_parms2$ivm_cov_par,
                                  ivm_min_age = ivm_parms2$ivm_min_age,
                                  ivm_max_age = ivm_parms2$ivm_max_age,
                                  IVRM_start = ivm_parms2$IVRM_start)

# Running each of the models wh0 (no ivermectin), wh1 (ivermectin, HR = 2 for 10 days) and
# wh2 (ivermectin, HR = 2 for 28 days)
res0 <- runfun(wh0)
res1 <- runfun(wh1)
res2 <- runfun(wh2)

# Plotting the results
cols <- c("grey40", "deeppink2", "deepskyblue3")
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))

# Clinical Incidence
plot(res1$t/365, res1$clin_inc0to80*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 2, col = cols[1], xlim = c(5.7, 10), ylim = c(0, 2000), las = 1)


res1$clin_inc0to80
res0$clin_inc0to80

res1$t
lines(res0$t/365, res0$clin_inc0to80*1000*365, lwd = 2, col = cols[2])
lines(res2$t/365, res2$clin_inc0to80*1000*365, lwd = 2, col = cols[3])
arrows(c(2250, 2280, 2310)/365, -50, c(2250, 2280, 2310)/365, 10, length = 0.15, lwd = 1, col = "goldenrod2")


# Slide Prevalence
plot(res0$t/365, res0$slide_prev0to80*100, type="l", ylab="Slide prevalence (%)", xlab="Year", lwd=3, col=cols[1],
     xlim = c(5.7, 8), ylim = c(0, 70), las=1)
lines(res1$t/365, res1$slide_prev0to80*100, lwd=3, col=cols[2])
lines(res2$t/365, res2$slide_prev0to80*100, lwd=3, col=cols[3])
arrows(c(2250, 2280, 2310)/365, -50, c(2250, 2280, 2310)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")

legend("topright", c("No endectocide", "IVM endectocide", "NTBC endectocide"),
       col = cols, lwd=3, bty="n", cex=0.8)

}


# Testing Plotting All Outputs
par(mfrow = c(1, 3))

# Clinical Incidence 0-5 years, 0 to 59 months, 0-80 years
plot(res0$t/365, res0$clin_inc0to5*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$clin_inc0to5*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$clin_inc0to5*1000*365, lwd = 3, col = cols[3])

plot(res0$t/365, res0$clin_inc0to59*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$clin_inc0to59*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$clin_inc0to59*1000*365, lwd = 3, col = cols[3])


# Slide Prevalence 0-5 years, 0 to 59 months, 0-80 years
plot(res0$t/365, res0$slide_prev0to5, type = "l", ylab = "Slide prevalence (%)", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1), las = 1)
lines(res1$t/365, res1$slide_prev0to5, lwd = 3, col = cols[2])
lines(res2$t/365, res2$slide_prev0to5, lwd = 3, col = cols[3])

plot(res0$t/365, res0$slide_prev0to59, type = "l", ylab = "Slide prevalence (%)", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1), las = 1)
lines(res1$t/365, res1$slide_prev0to59, lwd = 3, col = cols[2])
lines(res2$t/365, res2$slide_prev0to59, lwd = 3, col = cols[3])

plot(res0$t/365, res0$slide_prev0to80, type = "l", ylab = "Slide prevalence (%)", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1), las = 1)
lines(res1$t/365, res1$slide_prev0to80, lwd = 3, col = cols[2])
lines(res2$t/365, res2$slide_prev0to80, lwd = 3, col = cols[3])

# Deaths 0-5 years, 0 to 59 months, 0-80 years
plot(res0$t/365, res0$deaths_inc0to5*1000*365, type = "l", ylab = "Deaths", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$deaths_inc0to5*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$deaths_inc0to5*1000*365, lwd = 3, col = cols[3])

plot(res0$t/365, res0$deaths_inc0to59*1000*365, type = "l", ylab = "Deaths", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$deaths_inc0to59*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$deaths_inc0to59*1000*365, lwd = 3, col = cols[3])

plot(res0$t/365, res0$deaths_inc0to80*1000*365, type = "l", ylab = "Deaths", xlab = "Year",
     lwd = 3, col = cols[1], xlim = c(5.7, 8), ylim = c(0, 1000), las = 1)
lines(res1$t/365, res1$deaths_inc0to80*1000*365, lwd = 3, col = cols[2])
lines(res2$t/365, res2$deaths_inc0to80*1000*365, lwd = 3, col = cols[3])






# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# Create a vector of age categories for the model
init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 3 #low - 2, moderate - 15, high - 120

# Provide the length of time (in days) that you want to run the model for
time_period <- 3650 # run model for 10 years

# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

# Load ivermectin hazard
ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")

# Running the ivm_fun function to generate the extra endectocide specific parameters that you have to pass to the model
ivm_parms <- ivm_fun(IVM_start_times = 10000,            # time endectocide delivery occurs
                     time_period = time_period,         # time period for the model to run over
                     hazard_profile = ivm_haz$IVM_400_1_HS[1:23], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                     ivm_coverage = 0.8, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide

# Creates the odin model with all the required parameters - it is then ready to run
# Note this isn't the part where the model is actually run - that's below.
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

# This function runs the odin model with all the parameters that you previously generated with
# the ivRmectimn::create_r_model call
runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}

### The above code shows the workflow - the below code simulates the model and plots incidence and
### prevalence with and without ivermectin

#one scenario
ivm_parms0 <- ivm_fun(IVM_start_times = 10000, #no ivermectin: turning ivermectin on out of bounds of the model run time (3650 days)
                      time_period = time_period,
                      hazard_profile = ivm_haz$IVM_300_3_HS[1:23], #select 300 dosage here
                      ivm_coverage=0.8,
                      ivm_min_age=5,
                      ivm_max_age = 90)

wh0 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = 14,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms0$ttt,
                                  eff_len = ivm_parms0$eff_len,
                                  haz = ivm_parms0$haz,
                                  ivm_cov_par = ivm_parms0$ivm_cov_par,
                                  ivm_min_age = ivm_parms0$ivm_min_age,
                                  ivm_max_age = ivm_parms0$ivm_max_age,
                                  IVRM_start = ivm_parms0$IVRM_start)

#another scenario
ivm_parms1 <- ivm_fun(IVM_start_times = c(3120, 3150, 3180),
                      time_period = time_period,
                      hazard_profile = ivm_haz$IVM_300_3_HS[1:28],
                      ivm_coverage=0.8,
                      ivm_min_age=5,
                      ivm_max_age = 90)
wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = 14,
                                  country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                  admin2 = "Fatick",
                                  ttt = ivm_parms1$ttt,
                                  eff_len = ivm_parms1$eff_len,
                                  haz = ivm_parms1$haz,
                                  ivm_cov_par = ivm_parms1$ivm_cov_par,
                                  ivm_min_age = ivm_parms1$ivm_min_age,
                                  ivm_max_age = ivm_parms1$ivm_max_age,
                                  IVRM_start = ivm_parms1$IVRM_start)

# Running each of the models wh0 (no ivermectin), wh1 (ivermectin, HR = 2 for 10 days) and
# wh2 (ivermectin, HR = 2 for 28 days)
res0 <- runfun(wh0)
res1 <- runfun(wh1)

# Plotting the results
cols <- c("grey40", "deeppink2")
par(mfrow = c(1, 2), mar = c(5, 4, 1, 1))

# Clinical Incidence
plot(res0$t/365, res0$clin_inc0to80*1000*365, type = "l", ylab = "Annual incidence per 1,000", xlab = "Year",
     lwd = 2, col = cols[1], xlim = c(7.5, 10), ylim = c(0, 4000), las = 1)
lines(res1$t/365, res1$clin_inc0to80*1000*365, lwd = 2, col = cols[2])
arrows(c(3120, 3150, 3180)/365, -50, c(3120, 3150, 3180)/365, 10, length = 0.15, lwd = 1, col = "goldenrod2")


# Slide Prevalence
plot(res0$t/365, res0$slide_prev0to80*100, type="l", ylab="Slide prevalence (%)", xlab="Year", lwd=3, col=cols[1],
     xlim = c(7.5,10), ylim = c(0,70), las=1)
lines(res1$t/365, res1$slide_prev0to80*100, lwd=3, col=cols[2])
arrows(c(3120, 3150, 3180)/365, -50, c(3120, 3150, 3180)/365, 1, length = 0.1, lwd = 3, col = "goldenrod2")

legend("topright", c("No endectocide", "IVM endectocide"),
       col = cols, lwd=3, bty="n", cex=0.8)

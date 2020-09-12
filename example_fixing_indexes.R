# test new vector indexing
# to check HR=1 has no impact on malaria

rm(list = ls())


# Loading the ivRmectin package
devtools::load_all()

# Create a vector of age categories for the model
init_age <- c(0, 0.5, 1, 2, 3.5, 4, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# Provide a value of the annual EIR for this model run
init_EIR <- 10

# Provide the length of time (in days) that you want to run the model for
time_period <- 365 * 10 # run model for 7 years

# Sourcing the extra functions required to generate the endectocide specific parameters
source("R/mda_ivm_functions.R")

# Load ivermectin hazardd
ivm_haz <- read.table("data/ivermectin_hazards.txt", header=TRUE)

runfun <- function(mod_name){
  mod <- mod_name$generator(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period)
  op<- mod$transform_variables(modx)
  return(op)
}


# no ivermectin - ivm with HR = 1
ivm_parms_HR1 = ivm_fun(IVM_start_times = c(2190, 2250, 2310),
                     time_period = time_period,
                     hazard_profile = rep(1, 10),
                     ivm_coverage=0.8,
                     ivm_min_age=5,
                     ivm_max_age = 80)

# no ivermectin
ivm_parms_0 = ivm_fun(IVM_start_times = 10000,
                        time_period = time_period,
                        hazard_profile = rep(1, 10),
                        ivm_coverage=0.8,
                        ivm_min_age=5,
                        ivm_max_age = 80)

wh1 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide_fix_index.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  country = NULL,
                                  admin2 = NULL,
                                  ttt = ivm_parms_HR1$ttt,
                                  eff_len = ivm_parms_HR1$eff_len,
                                  haz = ivm_parms_HR1$haz,
                                  ivm_cov_par = ivm_parms_HR1$ivm_cov_par,
                                  ivm_min_age = ivm_parms_HR1$ivm_min_age,
                                  ivm_max_age = ivm_parms_HR1$ivm_max_age,
                                  IVRM_start = ivm_parms_HR1$IVRM_start)

wh2 <- ivRmectin:::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide_fix_index.R",
                                  num_int = 2,
                                  het_brackets = 5,
                                  age = init_age,
                                  init_EIR = init_EIR,
                                  country = NULL,
                                  admin2 = NULL,
                                  ttt = ivm_parms_0$ttt,
                                  eff_len = ivm_parms_0$eff_len,
                                  haz = ivm_parms_0$haz,
                                  ivm_cov_par = ivm_parms_0$ivm_cov_par,
                                  ivm_min_age = ivm_parms_0$ivm_min_age,
                                  ivm_max_age = ivm_parms_0$ivm_max_age,
                                  IVRM_start = ivm_parms_0$IVRM_start)

res1 <- runfun(wh1)
res2 <- runfun(wh2)

cols <- c("grey40", "deeppink2", "deepskyblue3")

par(mfrow = c(1,2))

plot(res1$t/365, res1$slide_prev0to80*100, type="l", ylab="Slide prevalence (%)",
     xlab="Year", lwd=3, col=cols[1], ylim = c(25, 30),
     xlim = c(5.7, 8), las=1, main = "Prevalence")
lines(res2$t/365, res2$slide_prev0to80*100, lwd=3, col=cols[2])

#abline(v = c(2190, 2250, 2310)/365)


plot(res1$t/365, res1$mv, type="l", ylab="Vector density",
     xlab="Year", lwd=3, col=cols[1], ylim = c(4, 6),
     xlim = c(5.7, 8), las=1, main = "Vector density")
lines(res2$t/365, res2$mv, lwd=3, col=cols[2])


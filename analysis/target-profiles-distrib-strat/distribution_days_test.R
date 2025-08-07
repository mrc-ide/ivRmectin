#understanding how distributions work

# Loading the ivRmectin package
devtools::load_all()
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)

time_period <- 365*3 # run model for 3y

# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

# Load ivermectin hazard
ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")

init_EIR <- 100

IVM_begin1 <- c(365)

IVM_mda_times <- c(IVM_begin1, IVM_begin1+10)

# Running the ivm_fun function to generate the extra endectocide specific parameters that you have to pass to the model
ivm_parms <- ivm_fun(IVM_start_times = IVM_mda_times,            # time endectocide delivery occurs
                     time_period = time_period,         # time period for the model to run over
                     hazard_profile = ivm_haz$IVM_300_3_HS[1:23], # dummy hazard profile - must be vector (we'll change this later on). for 400 dosage
                     ivm_coverage = 0.8, # proportion of population receiving the endectocide
                     ivm_min_age = 5, # youngest age group receiving endectocide
                     ivm_max_age = 90) # oldest age group receiving endectocide

# Creates the odin model with all the required parameters - it is then ready to run
# Note this isn't the part where the model is actually run - that's below.
wh <- ivRmectin::create_r_model(odin_model_path = "inst/extdata/odin_model_endectocide.R",
                                num_int = 1, # number of vector control (IRS and ITN) population groups
                                het_brackets = 5, # number of heterogeneous biting categories
                                #age = init_age, # the different age classes to be ran within the model
                                init_EIR = init_EIR, # the Entomological Innoculation Rate
                                #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
                                ttt = ivm_parms$ttt, # model specific parameter to control timing of endectocide delivery
                                eff_len = ivm_parms$eff_len, # number of days after receiving endectocide that HR is higher
                                haz = ivm_parms$haz, # hazard ratio for each off the eff_len number of days
                                ivm_cov_par = ivm_parms$ivm_cov_par, # proportion of popuulation receiving the endectocide
                                ivm_min_age = ivm_parms$ivm_min_age, # youngest age group receiving endectocide
                                ivm_max_age = ivm_parms$ivm_max_age, # oldest age group receiving endectocide
                                IVRM_start = ivm_parms$IVRM_start) # model specific parameter to control timing of endectocide delivery

out <- runfun(wh)
out_df <- as.data.frame(out)
ggplot(out_df, aes(x = t, y = mv))+
  geom_line()

ivm_parms$IVRM_start

ivm_haz_df <- as.data.frame(ivm_haz)
time <- seq(1, 90)

df_mda <- as.data.frame(time)
df_mda <- df_mda %>%
  mutate(coverage_MDA_1 = c(rep(0.4, 23),
                            rep(0, 90-23)),
                  coverage_MDA_2 = c(rep(0, 10),
                                     rep(0.3,23),
                                     rep(0, 90-23-10)),
         MDA_HAZ_1 = c(ivm_haz_df$IVM_300_3_HS[1:23], rep(0, 90-23)),
         MDA_HAZ_2 = c(rep(0, 10),
                       ivm_haz_df$IVM_300_3_HS[1:23],
                       rep(0, 90-23-10)),
         experienced_HAZ_1 = coverage_MDA_1*MDA_HAZ_1,
         experienced_HAZ_2 = coverage_MDA_2*MDA_HAZ_2)

ggplot(df_mda, aes(x = time, y = MDA_HAZ_1))+
  geom_line()+
  geom_line(data = df_mda, aes(x = time, y = MDA_HAZ_2), linetype = "dashed")+
  ylab("hazards in pop")

#40% of the population receive MDA 1 and 30% receive MDA 2



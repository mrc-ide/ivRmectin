#script to compare different model types
devtools::load_all()
require(tidyverse)
#how much does mu_h change across different values of Q0 and ivm_cov?

# Provide a value of the annual EIR for this model run
#init_EIR_vec <- c(2, 25, 100) #low - 2, moderate - 15, high - 120 --> Ellie: low = 2, med = 25, high = 100

# Provide the length of time (in days) that you want to run the model for
#time_period <- 3650 # run model for 10 years, turn ivermectin on when nets are 6m, 1y, 2.5yo
time_period <- 365*10
mda_int <- 30
ivm_cov_in = c(0.1, 0.9)

#ivm on when nets are 6 months old
#net_seq <- seq(100, 3650, by = 3*365) #if other interventions e.g. nets are on, will need to model IVM on time in relation to this
IVM_begin1 <- 365*2
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocidepecific parameters
source("R/mda_ivm_functions.R")

#itn_on <- 100 #introduce nets 100 days into simulation

net_seq <- seq(365, 3650, by = 3*365)
runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

eff_len2 <- 23
#here the hazard profile is 1 everyday
ivm_parms2 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, eff_len2),
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #redundant because not in the model
  ivm_min_age=5,
  ivm_max_age = 90)


mod1_run <-  ivRmectin:::create_r_model(
  odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
  num_int = 2,
  #num_int = 2,
  ITN_IRS_on = 365,
  itn_cov = 0.75,
  #het_brackets = 5,
  #age = init_age,
  init_EIR = 100,
  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
  #admin2 = "Fatick",
  ttt = ivm_parms1$ttt,
  eff_len = ivm_parms1$eff_len,
  haz = ivm_parms1$haz,
  ivm_cov_par = 0.8,
  ivm_min_age = ivm_parms1$ivm_min_age,
  ivm_max_age = ivm_parms1$ivm_max_age,
  IVRM_start = ivm_parms1$IVRM_start,
  Q0 = 0.9,
  init_ft = 0
)
df_mod1 <- as.data.frame(runfun(mod1_run)) %>%
  mutate(model = "hazards")

ggplot(df_mod1, aes(x = t, y = EIR_tot))+
  geom_line()

#fits to get endec_mu and wane
endec_mu_vec <- seq(0, 1, 0.1)
wane <-  seq(0, 0.1, 0.01)
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var2 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_2 <- list()
for (i in seq_len(nrow(df_var2))){
  my_list_2[[i]] <- as.numeric(df_var2[i,])
}

mod2_run <-  function(data_in){
  endec_mu_in <- data_in[1]
  wane_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2,
    ITN_IRS_on = 365,
    itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms2$ttt,
    eff_len = ivm_parms2$eff_len,
    haz = ivm_parms2$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms2$ivm_min_age,
    ivm_max_age = ivm_parms2$ivm_max_age,
    IVRM_start = ivm_parms2$IVRM_start,
    Q0 = 0.9,
    endec_mu = endec_mu_in,
    wane = wane_in,
    endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_mod2 <- function(){
  mod2_out_list <- lapply(my_list_2, mod2_run)
  res_mod2_out <- lapply(mod2_out_list, runfun)
  mod2_df <- do.call(rbind, sapply(1:(nrow(df_var2)), function(x){
    df <- as.data.frame(res_mod2_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, itn_cov))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "endec_mu"))}, simplify = F))
  return(mod2_df)
}

df_mod2 <- my_sim_mod2()

df_mod1_mda <- df_mod1 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
  select(t, mv, EIR_tot, ivm_cov, Ivtot)

df_mod2_mda <- df_mod2 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
  select(t, mv, EIR_tot, ref, Ivtot)

mod2_list <- split(df_mod2_mda, f = df_mod2_mda$ref)
error <- numeric()

#for (i in 1:length(endec_mu_vec)){
#  error2_low <- c(error2_low, sum((df_mod1_mda_low_cov$Ivtot - mod2_list[[i]]$Ivtot)^2))
#}

for (i in 1:nrow(df_var2)){
  error <- c(error, sum((df_mod1_mda$Ivtot - mod2_list[[i]]$Ivtot)^2))
}

range(error)
index <- which.min(error)
params <- df_var2[index,] #endec_mu = 0.2, wane = 0.03
index
out <- df_mod2 %>%
  filter(ref == index) %>%
  mutate(model = "malariasim-odin-exp-decay")

ggplot(out, aes(x = t, y = Ivtot))+
  geom_line()+
  xlim(365, 365*3)


ggplot(df_mod1, aes(x = t, y = Ivtot))+
  geom_line()


ggplot(out, aes(x = t, y = Ivtot))+
  geom_line()


df_mod1 <- df_mod1 %>%
  mutate(model = "odin-hazards")

ggplot(df_mod1, aes(x = t, y = Ivtot))+
  geom_line()+
  geom_line(data = out, aes(x = t, y = Ivtot), col = "blue")


#good overlap big big YAY!

ggplot(df_mod1, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_line(data = out, aes(x = t, y = EIR_tot), col = "blue")


ggplot(df_mod1, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  ylim(0.15, 0.8)+
  geom_line(data = out, aes(x = t, y = slide_prev0to5), col = "blue")




df_mod1_bind <- df_mod1 %>%
  select(t, mv, EIR_tot, slide_prev0to5, model, mu)

out_bind <- out %>%
  select(t, mv, EIR_tot, slide_prev0to5, model, mu)

models <- rbind(df_mod1_bind, out_bind)
write_rds(models, file = "analysis/exploring_interactions/malariasim-odin/hazards_malariasim_odin_nets.rds")


#fits to get endec_mu and wane
endec_mu_vec <- 0.2
wane <-  0.03
#endec_mu_vec <- 0.1
#wane <- 0.2
df_var2 <- expand.grid(endec_mu = endec_mu_vec, wane = wane)
#df_var2 <- expand.grid(endec_mu = 0.1)
my_list_2 <- list()
for (i in seq_len(nrow(df_var2))){
  my_list_2[[i]] <- as.numeric(df_var2[i,])
}



mod2_run2 <-  function(data_in){
  endec_mu_in <- data_in[1]
  wane_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2,
    ITN_IRS_on = 365,
    itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms2$ttt,
    eff_len = ivm_parms2$eff_len,
    haz = ivm_parms2$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms2$ivm_min_age,
    ivm_max_age = ivm_parms2$ivm_max_age,
    IVRM_start = ivm_parms2$IVRM_start,
    Q0 = 0.9,
    endec_mu = endec_mu_in,
    wane = wane_in,
    endec_on = IVM_start,
    init_ft = 0
  )
  return(output)
}


my_sim_mod2 <- function(){
  mod2_out_list <- lapply(my_list_2, mod2_run2)
  res_mod2_out <- lapply(mod2_out_list, runfun)
  mod2_df <- do.call(rbind, sapply(1:(nrow(df_var2)), function(x){
    df <- as.data.frame(res_mod2_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead, ivm_cov, Q0, endec_mu, EIR_tot, slide_prev0to5, Ivtot, wane, endec_killing, itn_cov))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model = "endec_mu"))}, simplify = F))
  return(mod2_df)
}

df_mod2 <- my_sim_mod2()

df_mod1_bind <- df_mod1 %>%
  select(t, mv, EIR_tot, slide_prev0to5, model, mu)

df_mod2_bind <- df_mod2 %>%
  select(t, mv, EIR_tot, slide_prev0to5, model, mu)

models <- rbind(df_mod1_bind, df_mod2_bind)

ggplot(models, aes(x = t, y = EIR_tot, col = as.factor(model)))+
  geom_line()

write_rds(models, file = "analysis/exploring_interactions/malariasim-odin/hazards_malariasim_odin_nets.rds")

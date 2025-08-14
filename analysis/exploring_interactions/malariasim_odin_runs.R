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
net_seq <- seq(100, 3650, by = 3*365)
IVM_begin1 <- net_seq[3]+180 # 6 months into new net distribution
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

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

#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

#here the hazard profile is 1 everyday
ivm_parms2 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start,
  time_period = time_period,
  hazard_profile = rep(1, 23),
  #hazard_profile = hazzy,
  ivm_coverage=ivm_cov,
  ivm_min_age=5,
  ivm_max_age = 90)

#get different combinations of ivm_cov and Q0
#ivm_cov_in <- seq(0.1, 0.9, by = 0.1)
#Q0_in <- seq(0.1, 0.9, by = 0.1)

ivm_cov_in <- c(0.1, 0.9)
Q0_in <- 0.9

df_var1 <- expand.grid(Q0 = Q0_in, ivm_cov = ivm_cov_in)
names(df_var1)

#df <- df[1,]
my_list <- list()
for (i in seq_len(nrow(df_var1))){
  my_list[[i]] <- as.numeric(df_var1[i,])
}

a <- Sys.time()

mod1 <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
  odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
  num_int = 1,
  #num_int = 2,
  #ITN_IRS_on = 100,
  #itn_cov = 0.75,
  #het_brackets = 5,
  #age = init_age,
  init_EIR = 100,
  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
  #admin2 = "Fatick",
  ttt = ivm_parms1$ttt,
  eff_len = ivm_parms1$eff_len,
  haz = ivm_parms1$haz,
  ivm_cov_par = ivm_cov_in,
  ivm_min_age = ivm_parms1$ivm_min_age,
  ivm_max_age = ivm_parms1$ivm_max_age,
  IVRM_start = ivm_parms1$IVRM_start,
  Q0 = Q0_in
)
  return(output)
}

my_sim_mod1 <- function(){
  mod1_out_list <- lapply(my_list, mod1)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var1)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, Q0, ivm_cov, slide_prev0to5, EIR_tot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "hazards"))}, simplify = F))
  return(mod1_df)
}

df_mod1 <- my_sim_mod1()

#fits to get endec_mu
endec_mu_vec <- seq(0, 0.5, 0.01)
df_var2 <- expand.grid(Q0 = Q0_in, ivm_cov= ivm_cov_in, endec_mu_in = endec_mu_vec)

my_list_2 <- list()
for (i in seq_len(nrow(df_var2))){
  my_list_2[[i]] <- as.numeric(df_var2[i,])
}
mod2 <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  endec_mu_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = 100,
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick",
    ttt = ivm_parms2$ttt,
    eff_len = ivm_parms2$eff_len,
    haz = ivm_parms2$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms2$ivm_min_age,
    ivm_max_age = ivm_parms2$ivm_max_age,
    IVRM_start = ivm_parms2$IVRM_start,
    Q0 = Q0_in,
    endec_mu = endec_mu_in
  )
  return(output)
}

my_sim_mod2 <- function(){
  mod2_out_list <- lapply(my_list_2, mod2)
  res_mod2_out <- lapply(mod2_out_list, runfun)
  mod2_df <- do.call(rbind, sapply(1:(nrow(df_var2)), function(x){
    df <- as.data.frame(res_mod2_out[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df, t, mu, mv, ivm_cov, Q0, mu_h))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, net_type = "mu_h"))}, simplify = F))
  return(mod2_df)
}

df_mod2 <- my_sim_mod2()

mod3 <- function(data_in){
  ivm_cov_in <- data_in[1]
  Q0_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
  odin_model_path = system.file("extdata/odin_model_endectocide_mu_h.R", package = "ivRmectin"),
  num_int = 1,
  #num_int = 2,
  #ITN_IRS_on = 100,
  #itn_cov = 0.75,
  #het_brackets = 5,
  #age = init_age,
  init_EIR = 100,
  #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
  #admin2 = "Fatick",
  ttt = ivm_parms2$ttt,
  eff_len = ivm_parms2$eff_len,
  haz = ivm_parms2$haz,
  ivm_cov_par = ivm_cov_in,
  ivm_min_age = ivm_parms2$ivm_min_age,
  ivm_max_age = ivm_parms2$ivm_max_age,
  IVRM_start = ivm_parms2$IVRM_start,
  Q0 = Q0_in,
  mu_h = 0.26
)
  return(output)
}

res_mod3 <- runfun(mod3)
df_mod3 <- as.data.frame(res_mod3)


df_mod1 <- df_mod1 %>%
  mutate(mod_type = "odin-hazards",
         endec_mu = "hazards") %>%
  select(t, mv, EIR_tot, slide_prev0to5, ivm_cov, Q0, mod_type)

df_mod2 <- df_mod2 %>%
  mutate(mod_type = "malariasim-odin") %>%
  select(t, mv, EIR_tot, slide_prev0to5, ivm_cov, Q0, mod_type)

df_mod3 <- df_mod3 %>%
  mutate(mod_type = "mu_h",
         endec_mu = "0.26") %>%
  select(t, mv, EIR_tot, slide_prev0to5, ivm_cov, Q0, mod_type)

mods <- do.call("rbind", list(df_mod1, df_mod2, df_mod3))

ggplot(mods, aes(x = t, y = mv, col = as.factor(mod_type)))+
  geom_line()

#fits to estimate mu_h or endec_mu
unique(df_mod2$ivm_cov)

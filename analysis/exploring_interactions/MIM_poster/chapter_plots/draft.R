#set up odin runs with odin_model_endectocide.R and exp decay model
time_period <- 365*10
mda_int <- 30

#introduce nets 1y into simulation

itn_on <- 365

#IVM MDA starts 1y after last campaign...should change this to 7 for 1y after
IVM_begin1 <- 365*6
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
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
ivm_parms2$eff_len
#get different combinations of ivm_cov and Q0
#ivm_cov_in <- seq(0.1, 0.9, by = 0.1)
#Q0_in <- seq(0.1, 0.9, by = 0.1)

#specify parameters




df_var1 <- readRDS("W:/endectocides-cluster/data/scenario-parameter-set.rds")
#one row
#df_var1 <- df_var1[1,]
my_list <- list()
for (i in seq_len(nrow(df_var1))){
  my_list[[i]] <- as.numeric(df_var1[i,])
}

a <- Sys.time()

mod1 <-  function(data_in){
  d_ITN0_in <- data_in[1]
  itn_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  bites_Bed_in <- data_in[4]
  ivm_cov_in <- data_in[5]
  Q0_in <- data_in[6]
  r_ITN0_in <- data_in[7]
  itn_half_life_in <- data_in[8]

  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 2,
    ITN_IRS_on = 365, #nets on 1y into sim
    itn_cov = itn_cov_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms1$ttt,
    eff_len = ivm_parms1$eff_len,
    haz = ivm_parms1$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms1$ivm_min_age,
    ivm_max_age = ivm_parms1$ivm_max_age,
    IVRM_start = ivm_parms1$IVRM_start

  )
  return(output)
}

my_sim_mod1 <- function(){
  mod1_out_list <- lapply(my_list, mod1)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var1)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIR_tot, Ivtot, IVRM_sr, bites_Bed))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "hazards"))}, simplify = F))
  return(mod1_df)
}

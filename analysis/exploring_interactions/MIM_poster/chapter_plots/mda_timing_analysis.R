devtools::load_all()
require(tidyverse)
require(ggpattern)
#set up odin runs with odin_model_endectocide.R and exp decay model
time_period <- 365*15
mda_int <- 30

#introduce nets 1y into simulation

itn_on <- 365*5

net_seq <- seq(itn_on, time_period, by = 3*365)


#IVM MDA starts 6m after last ITN campaign
IVM_begin1 <- net_seq[2]+180
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

#IVM MDA starts 1y after last ITN campaign
IVM_begin2 <- net_seq[2] + 365
IVM_start2 <- c(IVM_begin2, IVM_begin2+mda_int, IVM_begin2+mda_int+mda_int)

#IVM starts 2.5 after last ITN campaign
IVM_begin3 <- net_seq[2] + (2*365)
IVM_start3 <- c(IVM_begin3, IVM_begin3+mda_int, IVM_begin3+mda_int+mda_int)


ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
#itn_on <- 100 #introduce nets 100 days into simulation


runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#set IVM params for Hannah's mosq model with hazards#


ivm_parms1 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start1,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms2 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start2,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms3 <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = IVM_start3,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)


#specify parameters
df_var1 <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-mda-times.rds")

#one row
#df_var1 <- df_var1[1,]

my_list <- list()
for (i in seq_len(nrow(df_var1))){
  my_list[[i]] <- as.numeric(df_var1[i,])
}

a <- Sys.time()


#early MDA

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
    ITN_IRS_on = net_seq[1],
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
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "hazards-early-mda"))}, simplify = F))
  return(mod1_df)
}


df_var1_mda_early <- my_sim_mod1()
saveRDS(df_var1_mda_early, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_early.rds")


#medium MDA

mod2 <-  function(data_in){
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
    ITN_IRS_on = net_seq[1], #nets on 1y into sim
    itn_cov = itn_cov_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms2$ttt,
    eff_len = ivm_parms2$eff_len,
    haz = ivm_parms2$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms2$ivm_min_age,
    ivm_max_age = ivm_parms2$ivm_max_age,
    IVRM_start = ivm_parms2$IVRM_start

  )
  return(output)
}

my_sim_mod2 <- function(){
  mod1_out_list <- lapply(my_list, mod2)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var1)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "hazards-medium-mda"))}, simplify = F))
  return(mod1_df)
}


df_var1_mda_medium <- my_sim_mod2()
saveRDS(df_var1_mda_medium, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_medium.rds")

#late MDA

mod3 <-  function(data_in){
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
    ITN_IRS_on = net_seq[1], #nets on 1y into sim
    itn_cov = itn_cov_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms3$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start

  )
  return(output)
}

my_sim_mod3 <- function(){
  mod1_out_list <- lapply(my_list, mod3)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var1)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "hazards-late-mda"))}, simplify = F))
  return(mod1_df)
}


df_var1_mda_late <- my_sim_mod3()
saveRDS(df_var1_mda_late, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_late.rds")

#baseline scenario: no interventions
df_var_base <- df_var1 %>%
  select(init_EIR, Q0)

my_list_base <- list()
for (i in seq_len(nrow(df_var_base))){
  my_list_base[[i]] <- as.numeric(df_var_base[i,])
}


mod4 <-  function(data_in){
  init_EIR_in <- data_in[1]
  Q0_in <- data_in[2]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #ITN_IRS_on = net_seq[1], #nets on 1y into sim
    #itn_cov = itn_cov_in,
    #bites_Bed = bites_Bed_in,
    #d_ITN0 = d_ITN0_in,
    #r_ITN0 = r_ITN0_in,
    #itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms3$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start

  )
  return(output)
}

my_sim_mod4 <- function(){
  mod1_out_list <- lapply(my_list_base, mod4)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var_base)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "baseline"))}, simplify = F))
  return(mod1_df)
}


df_baseline <- my_sim_mod4()
saveRDS(df_baseline, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_baseline.rds")


df_baseline <- df_baseline %>%
  mutate(bites_Bed = 0.95, ivm_cov = 0.7) #to help binding


#nets only
mod5 <-  function(data_in){
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
    ITN_IRS_on = net_seq[1], #nets on 1y into sim
    itn_cov = itn_cov_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms3$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = 0, #no endectocide
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start

  )
  return(output)
}

my_sim_mod5 <- function(){
  mod1_out_list <- lapply(my_list, mod5)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var1)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "ITN-int"))}, simplify = F))
  return(mod1_df)
}

df_ITN <- my_sim_mod5()
saveRDS(df_ITN, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_ITN.rds")



df_timings <- do.call("rbind", list(df_baseline, df_ITN, df_var1_mda_early, df_var1_mda_medium, df_var1_mda_late))


#MDA timing analysis in odin_model_malariasim_exp_decay.R####

df_exp_decay_var <- df_var1 %>%
  select(-ivm_cov)

endec_mu_vec <- seq(0,1,0.1)
wane_vec <- seq(0,1,0.1)

malariasim_odin <- expand.grid(endec_mu_vec, wane_vec) %>%
  rename(endec_mu = Var1,
         wane = Var2)

df_exp_decay <- cross_join(df_exp_decay_var, malariasim_odin)

list_exp_decay <- list()

for (i in seq_len(nrow(df_exp_decay))){
  list_exp_decay[[i]] <- as.numeric(df_exp_decay[i,])
}

#early MDA

mod1b <-  function(data_in){
  d_ITN0_in <- data_in[1]
  itn_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  bites_Bed_in <- data_in[4]
  Q0_in <- data_in[5]
  r_ITN0_in <- data_in[6]
  itn_half_life_in <- data_in[7]
  endec_mu_in <- data_in[8]
  wane_in <- data_in[9]

  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    ITN_IRS_on = net_seq[1],
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
    ivm_cov_par = 0, #redundant parameter
    ivm_min_age = ivm_parms1$ivm_min_age,
    ivm_max_age = ivm_parms1$ivm_max_age,
    IVRM_start = ivm_parms1$IVRM_start,
    endec_mu = endec_mu_in,
    wane = wane_in

  )
  return(output)
}

my_sim_mod1b <- function(){
  mod1_out_list <- lapply(list_exp_decay, mod1b)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_exp_decay)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "malariasim-exp-early-mda"))}, simplify = F))
  return(mod1_df)
}


df_var1_mda_early_b <- my_sim_mod1b()
#saveRDS(df_var1_mda_early_b, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_early_b.rds")


#medium MDA

mod2b <-  function(data_in){
  d_ITN0_in <- data_in[1]
  itn_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  bites_Bed_in <- data_in[4]
  Q0_in <- data_in[5]
  r_ITN0_in <- data_in[6]
  itn_half_life_in <- data_in[7]
  endec_mu_in <- data_in[8]
  wane_in <- data_in[9]

  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    ITN_IRS_on = net_seq[1], #nets on 1y into sim
    itn_cov = itn_cov_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms2$ttt,
    eff_len = ivm_parms2$eff_len,
    haz = ivm_parms2$haz,
    ivm_cov_par = 0, #redundant parameter
    ivm_min_age = ivm_parms2$ivm_min_age,
    ivm_max_age = ivm_parms2$ivm_max_age,
    IVRM_start = ivm_parms2$IVRM_start,
    endec_mu = endec_mu_in,
    wane = wane_in

  )
  return(output)
}

my_sim_mod2b <- function(){
  mod1_out_list <- lapply(list_exp_decay, mod2b)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_exp_decay)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "malariasim-exp-medium-mda"))}, simplify = F))
  return(mod1_df)
}


df_var1_mda_medium_b <- my_sim_mod2b()
#saveRDS(df_var1_mda_medium_b, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_medium.rds")

#late MDA

mod3b <-  function(data_in){
  d_ITN0_in <- data_in[1]
  itn_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  bites_Bed_in <- data_in[4]
  Q0_in <- data_in[5]
  r_ITN0_in <- data_in[6]
  itn_half_life_in <- data_in[7]
  endec_mu_in <- data_in[8]
  wane_in <- data_in[9]

  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    ITN_IRS_on = net_seq[1], #nets on 1y into sim
    itn_cov = itn_cov_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms3$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = 0, #redundant parameter
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start,
    endec_mu = endec_mu_in,
    wane = wane_in

  )
  return(output)
}

my_sim_mod3b <- function(){
  mod1_out_list <- lapply(list_exp_decay, mod3b)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_exp_decay)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "malariasim-exp-late-mda"))}, simplify = F))
  return(mod1_df)
}


df_var1_mda_late_b <- my_sim_mod3b()
#saveRDS(df_var1_mda_late_b, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_late_b.rds")

#baseline scenario: no interventions
df_var_base <- df_var1 %>%
  select(init_EIR, Q0)

my_list_base <- list()
for (i in seq_len(nrow(df_var_base))){
  my_list_base[[i]] <- as.numeric(df_var_base[i,])
}


mod4b <-  function(data_in){
  init_EIR_in <- data_in[1]
  Q0_in <- data_in[2]
  endec_mu_in <- data_in[3]
  wane_in <- data_in[4]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 1,
    #ITN_IRS_on = net_seq[1], #nets on 1y into sim
    #itn_cov = itn_cov_in,
    #bites_Bed = bites_Bed_in,
    #d_ITN0 = d_ITN0_in,
    #r_ITN0 = r_ITN0_in,
    #itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms3$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = 0,
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start,
    endec_mu = endec_mu_in,
    wane = wane_in

  )
  return(output)
}

my_sim_mod4b <- function(){
  mod1_out_list <- lapply(my_list_base, mod4b)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_var_base)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "malariasim-exp-baseline"))}, simplify = F))
  return(mod1_df)
}


df_baselineb <- my_sim_mod4b()

df_baseline <- df_baseline %>%
  mutate(bites_Bed = 0.95, ivm_cov = 0.7) #to help binding


#nets only
mod5b <-  function(data_in){
  d_ITN0_in <- data_in[1]
  itn_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  bites_Bed_in <- data_in[4]
  Q0_in <- data_in[5]
  r_ITN0_in <- data_in[6]
  itn_half_life_in <- data_in[7]
  endec_mu_in <- data_in[8]
  wane_in <- data_in[9]

  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    ITN_IRS_on = net_seq[1], #nets on 1y into sim
    itn_cov = itn_cov_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in,
    Q0 = Q0_in,
    init_EIR = init_EIR_in,
    init_ft = 0,
    #ivm parameters
    ttt = ivm_parms3$ttt,
    eff_len = ivm_parms3$eff_len,
    haz = ivm_parms3$haz,
    ivm_cov_par = 0, #no endectocide
    ivm_min_age = ivm_parms3$ivm_min_age,
    ivm_max_age = ivm_parms3$ivm_max_age,
    IVRM_start = ivm_parms3$IVRM_start,
    endec_mu = endec_mu_in,
    wane = wane_in

  )
  return(output)
}

my_sim_mod5b <- function(){
  mod1_out_list <- lapply(list_exp_decay, mod5b)
  res_mod1_out <- lapply(mod1_out_list, runfun)
  mod1_df <- do.call(rbind, sapply(1:(nrow(df_exp_decay)), function(x){
    df <- as.data.frame(res_mod1_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, slide_prev0to5, EIRout, Ivtot, IVRM_sr, bites_Bed, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "malariasim-exp-decay-ITN-int"))}, simplify = F))
  return(mod1_df)
}

df_ITNb <- my_sim_mod5b()







df_timings %>%
  filter(t == net_seq[2]) %>%
  group_by(model_type, ref) %>%
  summarise(EIRout = EIRout) #ref 2 is the high EIR

pals <- c('#fbb4ae','#fed9a6','#ccebc5','#decbe4', '#b3cde3')

dynamics_timing_mda <- df_timings %>%
  filter(ref == 2) %>%
  mutate(model_type = factor(model_type,
                             levels = c("baseline", "ITN-int", "hazards-early-mda",
                                        "hazards-medium-mda", "hazards-late-mda"),
                             labels = c("Baseline", "ITN only", "Ivermectin MDA 6m after ITN campaign",
                                        "Ivermectin MDA 1y after ITN campaign", "Ivermectin MDA 2y after ITN campaign"))) %>%
  ggplot()+
  aes(x = (t-net_seq[2])/365, y = slide_prev0to5*100, col = model_type)+
  geom_line(linewidth = 2)+
  scale_colour_manual(values = pals, name = "Scenario")+
  theme_bw(base_size = 14)+
  coord_cartesian(ylim = c(0,100), xlim = c(0.15, 3))+
  xlab("Time (years) since last ITN campaign")+
  ylab("Slide prevalence (%) in children under 5-years-old")+
  theme(legend.position = c(0.6, 0.2))

#cases averted for each scenario
impact_baseline_early <- df_baseline %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_baseline = sum(clin_inc0to5))

impact_early_MDA <- df_var1_mda_early %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_early_MDA = sum(clin_inc0to5))

out_impact_early <- left_join(impact_baseline_early, impact_early_MDA) %>%
  mutate(case_avert_rel = ((tot_cases_baseline  - tot_cases_early_MDA)/tot_cases_baseline)*100)


#for medium
impact_baseline_medium <- df_baseline %>%
  filter(between(t, IVM_begin2, IVM_begin2+180)) %>% #choose same times, in case add on seasonal profiles
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_baseline = sum(clin_inc0to5))

impact_medium_MDA <- df_var1_mda_medium %>%
  filter(between(t, IVM_begin2, IVM_begin2+180)) %>%
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_medium_MDA = sum(clin_inc0to5))

out_impact_medium <- left_join(impact_baseline_medium, impact_medium_MDA) %>%
  mutate(case_avert_rel = ((tot_cases_baseline  - tot_cases_medium_MDA)/tot_cases_baseline)*100)

#late
impact_baseline_late <- df_baseline %>%
  filter(between(t, IVM_begin3, IVM_begin3+180)) %>% #choose same times, in case add on seasonal profiles
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_baseline = sum(clin_inc0to5))

impact_late_MDA <- df_var1_mda_late %>%
  filter(between(t, IVM_begin3, IVM_begin3+180)) %>%
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_late_MDA = sum(clin_inc0to5))

out_impact_late <- left_join(impact_baseline_late, impact_late_MDA) %>%
  mutate(case_avert_rel= ((tot_cases_baseline  - tot_cases_late_MDA)/tot_cases_baseline)*100)


#suggest combined endec + ITN early in campaign is better (relative to baseline) than later
#this is picking up the net effect; the net itself is better earlier in the campaign then later

#to isolate impact of the endectocide only, compare endec+ITN to ITN only scenario
impact_ITN_early <- df_ITN %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_ITN= sum(clin_inc0to5))

out_endec_impact_early <- left_join(impact_ITN_early,impact_early_MDA) %>%
  mutate(case_avert_rel = ((tot_cases_ITN  - tot_cases_early_MDA)/tot_cases_ITN)*100)


impact_ITN_medium <- df_ITN %>%
  filter(between(t, IVM_begin2, IVM_begin2+180)) %>%
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_ITN = sum(clin_inc0to5))

out_endec_impact_medium <- left_join(impact_ITN_medium,impact_medium_MDA) %>%
  mutate(case_avert_rel = ((tot_cases_ITN  - tot_cases_medium_MDA)/tot_cases_ITN)*100)

impact_ITN_late <- df_ITN %>%
  filter(between(t, IVM_begin3, IVM_begin3+180)) %>%
  group_by(ref) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  summarise(tot_cases_ITN = sum(clin_inc0to5))

out_endec_impact_late <- left_join(impact_ITN_late,impact_late_MDA) %>%
  mutate(case_avert_rel = ((tot_cases_ITN  - tot_cases_late_MDA)/tot_cases_ITN)*100)

#compared to ITN only scenarios, the impact of the endectocide is consistent regardless of when the distribution is done

#compared to baseline scenario, impact of interventions is greater when MDA is early (likely because the impact of the net is higher here, so bigger net effect)

output_impact_early <- out_impact_early %>%
  select(ref, case_avert_rel) %>%
  mutate(scenario = "early-MDA",
         measure = "rel-baseline")
output_impact_medium <- out_impact_medium %>%
  select(ref, case_avert_rel) %>%
  mutate(scenario = "medium-MDA",
         measure = "rel-baseline")
output_impact_late <- out_impact_late %>%
  select(ref, case_avert_rel) %>%
  mutate(scenario = "late-MDA",
         measure = "rel-baseline")

output_impact_early_MDA <- out_endec_impact_early %>%
  select(ref, case_avert_rel) %>%
  mutate(scenario = "early-MDA",
         measure = "rel-ITN")
output_impact_medium_MDA <- out_endec_impact_medium %>%
  select(ref, case_avert_rel) %>%
  mutate(scenario = "medium-MDA",
         measure = "rel-ITN")
output_impact_late_MDA <- out_endec_impact_late %>%
  select(ref, case_avert_rel) %>%
  mutate(scenario = "late-MDA",
         measure = "rel-ITN")

output_impact_rel_baseline <- do.call("rbind", list(output_impact_early,output_impact_medium,output_impact_late,
                                                    output_impact_early_MDA, output_impact_medium_MDA, output_impact_late_MDA))



pals2 <- c('#ccebc5','#decbe4', '#b3cde3')

case_avert_plot <- output_impact_rel_baseline %>%
  mutate(scenario = factor(scenario,
                             levels = c("early-MDA", "medium-MDA", "late-MDA"),
                             labels = c("6m after ITN campaign", "1y after MDA campaign", "2y after MDA campaign"))) %>%
  filter(ref == 2) %>%
  ggplot()+
  aes(x = scenario, y = case_avert_rel, fill = scenario, pattern = measure)+
  geom_bar_pattern(stat = "identity",
                   position = position_dodge(),
                   colour = "black",                # Border of bars
                   pattern_colour = "black",        # Pattern line color
                   pattern_fill = NA,               # Transparent so bar fill shows
                   pattern_density = 0.4,
                   pattern_spacing = 0.05,
                   pattern_key_scale_factor = 0.5)+
  scale_pattern_manual(values = c("none", "stripe"), name = "Efficacy measurement",
                       labels = c("Relative to baseline (no interventions) scenario",
                                  "Relative to scenario with historic use of ITNs"))+
  theme_bw(base_size = 14)+
  coord_cartesian(ylim = c(0,100))+
  xlab("Timing of MDA in relation to ITN campaign")+
  ylab("Cases averted (%) in children under 5-years-old")+
  scale_fill_manual(values = pals2)+
  theme(legend.position = c(0.6, 0.8))+
  guides(fill = "none", pattern_spacing = 0.5,
         pattern = guide_legend(
           override.aes = list(fill = "white")))

mda_timing_plot_fig2 <- cowplot::plot_grid(dynamics_timing_mda, case_avert_plot,
                   labels = c("A", "B"))
ggsave(mda_timing_plot, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/mda_timing_plot_fig2.pdf")


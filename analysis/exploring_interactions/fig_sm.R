#supp figure to look at

bb_antag_itn_ivm <- readRDS("X:/endectocides-cluster/data/antag_LLIN_IVM_bites_Bed_loop.rds")
bb_antag_itn <- readRDS("X:/endectocides-cluster/data/antag_LLIN_bites_Bed_loop.rds")

bb_add_itn_ivm <- readRDS("X:/endectocides-cluster/data/add_LLIN_IVM_bites_Bed_loop.rds")
bb_add_itn <- readRDS("X:/endectocides-cluster/data/add_LLIN_bites_Bed_loop.rds")

#extract baseline from antag

EIR_vals <- bb_antag_itn_ivm %>%
  group_by(bites_Bed, itn_cov, d_ITN0) %>% #changed species, itn coverage and resistance
  filter(t == 1)
EIR_vals2 <- unique(EIR_vals$EIR_tot) #2.781299  34.766235 139.064942. use this to ref with the EIRs

fillNAgaps <- function(x, firstBack=FALSE) {
  ## NA's in a vector or factor are replaced with last non-NA values
  ## If firstBack is TRUE, it will fill in leading NA's with the first
  ## non-NA value. If FALSE, it will not change leading NA's.

  # If it's a factor, store the level labels and convert to integer
  lvls <- NULL
  if (is.factor(x)) {
    lvls <- levels(x)
    x    <- as.integer(x)
  }

  goodIdx <- !is.na(x)

  # These are the non-NA values from x only
  # Add a leading NA or take the first good value, depending on firstBack
  if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
  else             goodVals <- c(NA,            x[goodIdx])

  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1

  x <- goodVals[fillIdx]

  # If it was originally a factor, convert it back
  if (!is.null(lvls)) {
    x <- factor(x, levels=seq_along(lvls), labels=lvls)
  }

  x
}

bb_antag_itn_ivm <- bb_antag_itn_ivm %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))

bb_antag_itn_ivm$init_EIR <- fillNAgaps(bb_antag_itn_ivm$init_EIR)

bb_antag_itn <- bb_antag_itn %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))

bb_antag_itn$init_EIR <- fillNAgaps(bb_antag_itn$init_EIR)

bb_add_itn_ivm <- bb_add_itn_ivm %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))
bb_add_itn_ivm$init_EIR <- fillNAgaps(bb_add_itn_ivm$init_EIR)

bb_add_itn <- bb_add_itn %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))
bb_add_itn$init_EIR <- fillNAgaps(bb_add_itn$init_EIR)

#check
unique(bb_antag_itn_ivm$init_EIR)
unique(bb_antag_itn$init_EIR)
unique(bb_add_itn_ivm$init_EIR)
unique(bb_add_itn$init_EIR)

#plotting the mean prevalence and EIR in the ivermectin distribution period
IVM_begin <- (365*8)+180
mda_int <- 30

IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)

eff_len <- 23
ivm_on <- IVM_start[1] #3100
ivm_off <- IVM_start[3]+eff_len #3183


#df_antag <- rbind(df_antag_ITN, df_antag_ITN_IVM)
#df_add <- rbind(df_add_ITN, df_add_ITN_IVM)

#extract the baseline values from each model type (from _ITN only)
bb_baseline_vals_antag <- bb_antag_itn %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  group_by(bites_Bed, itn_cov, d_ITN0, init_EIR) %>%
  summarise(mean_prev_baseline = mean(slide_prev0to5),
            mean_EIR_baseline = mean(EIR_tot),
            mean_avhc_baseline = mean(avhc),
            mean_mu_baseline = mean(mu))

bb_antag_eff <- bb_antag_itn_ivm %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  group_by(bites_Bed, itn_cov, d_ITN0, init_EIR) %>%
  summarise(mean_prev = mean(slide_prev0to5),
            mean_EIR = mean(EIR_tot),
            mean_avhc = mean(avhc),
            mean_mu = mean(mu)) %>%
  left_join(bb_baseline_vals_antag) %>%
  mutate(eff_prev_antag = mean_prev_baseline - mean_prev,
         eff_eir_antag = mean_EIR_baseline - mean_EIR) %>%
  select(species, itn_cov, d_ITN0, init_EIR, eff_prev_antag, eff_eir_antag)
#repeat this for the additive model
baseline_vals_add <- df_add_ITN %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  group_by(species, itn_cov, d_ITN0, init_EIR) %>%
  summarise(mean_prev_baseline = mean(slide_prev0to5),
            mean_EIR_baseline = mean(EIR_tot),
            mean_avhc_baseline = mean(avhc),
            mean_mu_baseline = mean(mu))

df_add_eff <- df_add_ITN_IVM %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  group_by(species, itn_cov, d_ITN0, init_EIR) %>%
  summarise(mean_prev = mean(slide_prev0to5),
            mean_EIR = mean(EIR_tot),
            mean_avhc = mean(avhc),
            mean_mu = mean(mu)) %>%
  left_join(baseline_vals_antag) %>%
  mutate(eff_prev_add = mean_prev_baseline - mean_prev,
         eff_eir_add = mean_EIR_baseline - mean_EIR) %>%
  select(species, itn_cov, d_ITN0, init_EIR, eff_prev_add, eff_eir_add)

#rbind them
df_eff <- left_join(df_antag_eff, df_add_eff)
df_eff$init_EIR <- factor(df_eff$init_EIR, levels = c("low", "medium", "high"))


require(tidyverse)
#comparing the odin and odin exp decay output

IVM_begin1 <- 365*6
mda_int <- 30
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)


df_mod1 <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/output_mod1.rds")
df_mod2 <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/odin_exp_decay_compare.rds")

#save locally
saveRDS(df_mod1, file = "analysis/exploring_interactions/malariasim-odin/output_mod1.rds")
saveRDS(df_mod2, file = "analysis/exploring_interactions/malariasim-odin/odin_exp_decay_compare.rds")

#read in the parameter sets to join on any extra information
param_set_df_mod1 <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set.rds")
param_set_df_mod2 <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-exp-decay.rds")


param_set_df_mod1 <- param_set_df_mod1 %>%
  mutate(ref = row_number())

param_set_df_mod2 <- param_set_df_mod2 %>%
  mutate(ref = row_number())

covs <- unique(df_mod1$ivm_cov)

df_mod1 <- df_mod1 %>%
  mutate(ivm_cov = case_when(ivm_cov == covs[1] ~ 0.5,
                             ivm_cov == covs[2] ~ 0.9))


df_mod1 <- left_join(df_mod1, param_set_df_mod1, by = "ref")
df_mod2 <- left_join(df_mod2, param_set_df_mod2, by = "ref")


df_mod1_mda <- df_mod1 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
  select(t, mv, init_EIR, EIR_tot, bites_Bed.x, ivm_cov.x, Ivtot, ref) %>%
  rename(ivm_cov = ivm_cov.x,
         bites_Bed = bites_Bed.x)

df_mod2_mda <- df_mod2 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
  select(t, mv, init_EIR, EIR_tot, bites_Bed.x, Ivtot, ref, endec_mu.x, wane.x) %>%
  rename(endec_mu = endec_mu.x,
         wane_endec = wane.x,
         bites_Bed = bites_Bed.x)

library(dplyr)

# Define unique values of `init_EIR` and `bites_Bed`
EIR_in <- unique(df_mod1_mda$init_EIR)
bites_Bed_in <- unique(df_mod1_mda$bites_Bed)

# Define parameter grids
endec_mu_vec <- seq(0, 1, 0.1)
wane_vec <- seq(0, 0.1, 0.01)
combos <- expand.grid(endec_mu = endec_mu_vec, wane_endec = wane_vec)

# Split data by coverage
data_splits <- list(
  "high_cov" = df_mod1_mda %>% filter(ivm_cov == 0.9),
  "low_cov" = df_mod1_mda %>% filter(ivm_cov == 0.5)
)




# Initialize results storage
results <- list()

# Iterate over data splits, `init_EIR`, and `bites_Bed`
for (split_name in names(data_splits)) {
  df_split <- data_splits[[split_name]]

  for (eir_val in EIR_in) {
    for (bites_val in bites_Bed_in) {

      # Filter the current split for specific `init_EIR` and `bites_Bed`
      df_1 <- df_split %>%
        filter(init_EIR == eir_val, bites_Bed == bites_val)

      # Skip if no data for this combination
      if (nrow(df_1) == 0) next

      # Filter df_mod2_mda for the same combination
      df_2 <- df_mod2_mda %>%
        filter(init_EIR == eir_val, bites_Bed == bites_val)

      # Skip if no data in df_mod2_mda
      if (nrow(df_2) == 0) next

      # Split df_2 by `ref`
      mod2_list <- split(df_2, f = df_2$ref)

      # Compute error for each combo of endec_mu and wane_endec
      error <- numeric()
      for (i in 1:nrow(combos)) {
        error <- c(error, sum((df_1$Ivtot - mod2_list[[i]]$Ivtot)^2))
      }

      # Find the best combination
      index <- which.min(error)
      best_fit <- combos[index, ]

      # Store results for this combination of `init_EIR` and `bites_Bed`
      results[[paste(split_name, "EIR", eir_val, "bitesBed", bites_val, sep = "_")]] <- list(
        split = split_name,
        init_EIR = eir_val,
        bites_Bed = bites_val,
        best_fit = best_fit,
        error = error[index]
      )
    }
  }
}

# Combine results into a dataframe for easier analysis
results_df <- do.call(rbind, lapply(names(results), function(name) {
  data.frame(
    scenario = name,
    split = results[[name]]$split,
    init_EIR = results[[name]]$init_EIR,
    bites_Bed = results[[name]]$bites_Bed,
    endec_mu = results[[name]]$best_fit$endec_mu,
    wane_endec = results[[name]]$best_fit$wane_endec,
    error = results[[name]]$error
  )
}))

# View the results
print(results_df)

######

mod2_high_covs_EIR_30_Bites_0.481 <- df_mod2 %>%
  filter(bites_Bed.x == bites_Bed_in[1], init_EIR == EIR_in[1], endec_mu.x == 0.4, wane.x == 0.06) %>%
  select(-c(endec_mu.x, wane.x)) %>%
  select(t, EIR_tot, slide_prev0to5, Ivtot, init_EIR, bites_Bed.x) %>%
  mutate(model = "endec_mu_decay") %>%
  rename(bites_Bed = bites_Bed.x)
names(mod2_high_covs_EIR_30_Bites_0.481)

mod1_high_covs_EIR_30_Bites_0.481 <- df_mod1 %>%
  filter(bites_Bed.x == bites_Bed_in[1], init_EIR == EIR_in[1], ivm_cov.y == 0.9) %>%
  select(-c(ivm_cov.y)) %>%
  select(t, EIR_tot, slide_prev0to5, Ivtot, init_EIR, bites_Bed.x) %>%
  rename(bites_Bed = bites_Bed.x) %>%
  mutate(model = "hazards")
names(mod1_high_covs_EIR_30_Bites_0.481)

mod2_high_covs_EIR_30_Bites_0.988 <- df_mod2 %>%
  filter(bites_Bed.x == bites_Bed_in[2], init_EIR == EIR_in[1], endec_mu.x == 0.4, wane.x == 0.06) %>%
  select(-c(endec_mu.x, wane.x)) %>%
  select(t, EIR_tot, slide_prev0to5, Ivtot, init_EIR, bites_Bed.x) %>%
  mutate(model = "endec_mu_decay") %>%
  rename(bites_Bed = bites_Bed.x)
names(mod2_high_covs_EIR_30_Bites_0.988)

mod1_high_covs_EIR_30_Bites_0.988 <- df_mod1 %>%
  filter(bites_Bed.x == bites_Bed_in[2], init_EIR == EIR_in[1], ivm_cov.y == 0.9) %>%
  select(-c(ivm_cov.y)) %>%
  select(t, EIR_tot, slide_prev0to5, Ivtot, init_EIR, bites_Bed.x) %>%
  rename(bites_Bed = bites_Bed.x) %>%
  mutate(model = "hazards")
names(mod1_high_covs_EIR_30_Bites_0.988)

combo_list <- list(mod2_high_covs_EIR_30_Bites_0.481,
                   mod1_high_covs_EIR_30_Bites_0.481,
                   mod2_high_covs_EIR_30_Bites_0.988,
                   mod1_high_covs_EIR_30_Bites_0.988
                   )

combo_high_covs_EIR_30 <- do.call("rbind", combo_list)

ggplot(combo_high_covs_EIR_30, aes(x = t, y = Ivtot, col = as.factor(model)))+
  geom_line()+
  facet_grid(~bites_Bed, labeller = label_both)



combo_high_covs_EIR_30_Bites_0.481 <- rbind(mod1_high_covs_EIR_30_Bites_0.481, mod2_high_covs_EIR_30_Bites_0.481)

ggplot(combo_high_covs_EIR_30_Bites_0.481, aes(x = t, y = Ivtot, col = as.factor(model)))+
  geom_line()

ggplot(combo_high_covs_EIR_30_Bites_0.481, aes(x = t, y = EIR_tot, col = as.factor(model)))+
  geom_line()
ggplot(combo_high_covs_EIR_30_Bites_0.481, aes(x = t, y = slide_prev0to5, col = as.factor(model)))+
  geom_line()











#for low coverage
mod2_low_covs_EIR_30_Bites_0.481 <- df_mod2 %>%
  filter(bites_Bed.x == bites_Bed_in[1], init_EIR == EIR_in[1], endec_mu.x == 0.2, wane.x == 0.02) %>%
  select(-c(endec_mu.x, wane.x)) %>%
  select(t, EIR_tot, slide_prev0to5, Ivtot, init_EIR, bites_Bed.x) %>%
  mutate(model = "endec_mu_decay") %>%
  rename(bites_Bed = bites_Bed.x)
names(mod2_low_covs_EIR_30_Bites_0.481)

mod1_low_covs_EIR_30_Bites_0.481 <- df_mod1 %>%
  filter(bites_Bed.x == bites_Bed_in[1], init_EIR == EIR_in[1], ivm_cov.y == 0.5) %>%
  select(-c(ivm_cov.y)) %>%
  select(t, EIR_tot, slide_prev0to5, Ivtot, init_EIR, bites_Bed.x) %>%
  rename(bites_Bed = bites_Bed.x) %>%
  mutate(model = "hazards")
names(mod1_low_covs_EIR_30_Bites_0.481)

combo_low_covs_EIR_30_Bites_0.481 <- rbind(mod1_low_covs_EIR_30_Bites_0.481, mod2_low_covs_EIR_30_Bites_0.481)

ggplot(combo_low_covs_EIR_30_Bites_0.481, aes(x = t, y = Ivtot, col = as.factor(model)))+
  geom_line()

ggplot(combo_low_covs_EIR_30_Bites_0.481, aes(x = t, y = EIR_tot, col = as.factor(model)))+
  geom_line()

ggplot(combo_low_covs_EIR_30_Bites_0.481, aes(x = t, y = slide_prev0to5, col = as.factor(model)))+
  geom_line()

#ivm_cov_values <- unique(df_mod1$ivm_cov)
#
#data_splits <- list()
#for (cov in ivm_cov_values) {
#  #filter by coverage
#  filtered_data <- df_mod1_mda %>%
#    filter(ivm_cov == cov)
#
#  #generate names for the list element, based on coverage
#  cov_name <- paste0("cov_", cov)
#
#  #add the filtered data to the list
#  data_splits[[cov_name]] <- filtered_data
#}
##then print the list names and first few rows for each filtered dataframe
#for (cov_name in names(data_splits)){
#  print(paste("Data for", cov_name))
#  print(head(data_splits[[cov_name]]))
#}



##manual check for that mortality rate.
EIR_in <- unique(df_mod1_mda$init_EIR)
bites_Bed_in <- unique(df_mod1_mda$bites_Bed)

endec_mu_vec <- seq(0, 1, 0.1)
wane_vec <- seq(0,0.1,0.01)

combos <- expand.grid(endec_mu = endec_mu_vec, wane_endec = wane_vec)


#split by coverage first
df_mod1_mda_low_cov <- df_mod1_mda %>%
  filter(ivm_cov == 0.5)

df_mod1_mda_high_cov <- df_mod1_mda %>%
  filter(ivm_cov == 0.9)

df_1_A <- df_mod1_mda_high_cov %>%
  filter(init_EIR == EIR_in[1],
         bites_Bed == bites_Bed_in[2])

ggplot(df_1_A, aes(x = t, y = Ivtot))+
  geom_line()

df_2_A <- df_mod2_mda %>%
  filter(init_EIR == EIR_in[1],
         bites_Bed == bites_Bed_in[2])

ggplot(df_2_A, aes(x = t, y = Ivtot))+
  geom_line()

mod2_list <- split(df_2_A, f = df_2_A$ref) #split for each wane and endec_mu
errorA <- numeric()
for(i in 1:nrow(combos)){
  errorA <- c(errorA, sum((df_1_A$Ivtot - mod2_list[[i]]$Ivtot)^2))
}
index_A <- which.min(errorA)
result <- combos[index_A,] #endec_mu = 0.2, wane_endec = 0.1



#check no intervention
df_mod3 <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/output_mod1_no_int.rds")
df_mod4 <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/odin_exp_decay_compare_no_int.rds")

saveRDS(df_mod3, file = "analysis/exploring_interactions/malariasim-odin/output_mod1_no_int.rds")
saveRDS(df_mod4, file = "analysis/exploring_interactions/malariasim-odin/odin_exp_decay_compare_no_int.rds")


param_set_df_mod3 <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-no-int.rds")
param_set_df_mod4 <- readRDS("analysis/exploring_interactions/malariasim-odin/scenario-parameter-set-exp-decay-no-int.rds")


df_mod3 <- df_mod3 %>%
  mutate(ivm_cov = case_when(ivm_cov == covs[1] ~ 0.5,
                             ivm_cov == covs[2] ~ 0.9))

param_set_df_mod3 <- param_set_df_mod3 %>%
  mutate(ref = row_number())

param_set_df_mod4 <- param_set_df_mod4 %>%
  mutate(ref = row_number())

df_mod3 <- left_join(df_mod3, param_set_df_mod3, by = "ref")
df_mod4 <- left_join(df_mod4, param_set_df_mod4, by = "ref")


covs <- unique(df_mod3$ivm_cov.x)

df_mod3_mda <- df_mod3 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
  select(t, mv, init_EIR, EIR_tot, ivm_cov.x, Ivtot, ref) %>%
  rename(ivm_cov = ivm_cov.x)

df_mod4_mda <- df_mod4 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
  select(t, mv, init_EIR, EIR_tot, Ivtot, endec_mu.x, wane.x, ref) %>%
  rename(endec_mu = endec_mu.x,
         wane_endec = wane.x)


# Define parameter grids
endec_mu_vec <- seq(0, 1, 0.1)
wane_vec <- seq(0, 0.1, 0.01)
combos <- expand.grid(endec_mu = endec_mu_vec, wane_endec = wane_vec)

# Split data by coverage
data_splits <- list(
  "high_cov" = df_mod3_mda %>% filter(ivm_cov == 0.9),
  "low_cov" = df_mod3_mda %>% filter(ivm_cov == 0.5)
)

# Initialize results storage
results_no_int <- list()

# Iterate over data splits, `init_EIR`, and `bites_Bed`
for (split_name in names(data_splits)) {
  df_split <- data_splits[[split_name]]

  for (eir_val in EIR_in) {

      # Filter the current split for specific `init_EIR`
      df_1 <- df_split %>%
        filter(init_EIR == eir_val)

      # Skip if no data for this combination
      if (nrow(df_1) == 0) next

      # Filter df_mod4_mda for the same combination
      df_4 <- df_mod4_mda %>%
        filter(init_EIR == eir_val)

      # Skip if no data in df_mod2_mda
      if (nrow(df_4) == 0) next

      # Split df_2 by `ref`
      mod2_list <- split(df_4, f = df_4$ref)

      # Compute error for each combo of endec_mu and wane_endec
      error <- numeric()
      for (i in 1:nrow(combos)) {
        error <- c(error, sum((df_1$Ivtot - mod2_list[[i]]$Ivtot)^2))
      }

      # Find the best combination
      index <- which.min(error)
      best_fit <- combos[index, ]

      # Store results for this combination of `init_EIR`
      results_no_int[[paste(split_name, "EIR", eir_val, sep = "_")]] <- list(
        split = split_name,
        init_EIR = eir_val,
        best_fit = best_fit,
        error = error[index]
      )

  }
}

# Combine results into a dataframe for easier analysis
results_df_no_int <- do.call(rbind, lapply(names(results_no_int), function(name) {
  data.frame(
    scenario = name,
    split = results_no_int[[name]]$split,
    init_EIR = results_no_int[[name]]$init_EIR,
    endec_mu = results_no_int[[name]]$best_fit$endec_mu,
    wane_endec = results_no_int[[name]]$best_fit$wane_endec,
    error = results_no_int[[name]]$error
  )
}))

# View the results
print(results_df_no_int)


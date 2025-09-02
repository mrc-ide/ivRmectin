require(tidyverse)
#comparing the odin and odin exp decay output

IVM_begin1 <- 365*6
mda_int <- 30
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

#from model_runs_odin_exp_decay.R


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
  select(t, mv, init_EIR, EIR_tot, bites_Bed.x, ivm_cov.x, Ivtot, d_ITN0, r_ITN0, itn_half_life, itn_cov, Q0.x, ref) %>%
  rename(ivm_cov = ivm_cov.x,
         bites_Bed = bites_Bed.x,
         Q0 = Q0.x)

df_mod2_mda <- df_mod2 %>%
  filter(between(t, IVM_start[1], IVM_start[3]+23)) %>%
  select(t, mv, init_EIR, EIR_tot, bites_Bed.x, Ivtot, ref, endec_mu.x, wane.x, d_ITN0, r_ITN0, itn_half_life, itn_cov, Q0.x) %>%
  rename(endec_mu = endec_mu.x,
         wane_endec = wane.x,
         bites_Bed = bites_Bed.x,
         Q0 = Q0.x)



# Define unique values of `init_EIR` and `bites_Bed` and resistance levels
EIR_in <- unique(df_mod1_mda$init_EIR)
bites_Bed_in <- unique(df_mod1_mda$bites_Bed)
Q0_in <- unique(df_mod1_mda$Q0)
d_ITN0_in <- unique(df_mod1_mda$d_ITN0)
itn_cov_in <- unique(df_mod1_mda$itn_cov)

# Define parameter grids
endec_mu_vec <- seq(0.01, 1, 0.11)
wane_vec <- seq(0.24,0.70, length.out = 5)
combos <- expand.grid(endec_mu = endec_mu_vec, wane_endec = wane_vec)

# Split data by coverage
data_splits <- list(
  "high_cov" = df_mod1_mda %>% filter(ivm_cov == 0.9),
  "low_cov" = df_mod1_mda %>% filter(ivm_cov == 0.5)
)

# Initialize results storage
results <- list()

# Iterate over data splits, `init_EIR`, and `bites_Bed`, res, itn_cov and Q0
for (split_name in names(data_splits)) {
  df_split <- data_splits[[split_name]]

  for (eir_val in EIR_in) {
    for (bites_val in bites_Bed_in) {
      for(res_val in d_ITN0_in){
        for(cov_val in itn_cov_in){
          for(Q0_val in Q0_in){

      # Filter the current split for specific `init_EIR` and `bites_Bed` and resistance
      df_1 <- df_split %>%
        filter(init_EIR == eir_val, bites_Bed == bites_val, d_ITN0 == res_val, itn_cov == cov_val, Q0 == Q0_val)

      # Skip if no data for this combination
      if (nrow(df_1) == 0) next

      # Filter df_mod2_mda for the same combination
      df_2 <- df_mod2_mda %>%
        filter(init_EIR == eir_val, bites_Bed == bites_val, d_ITN0 == res_val, itn_cov == cov_val, Q0 == Q0_val)

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

      # Store results for this combination of `init_EIR` and `bites_Bed` and resistance
      results[[paste(split_name, "EIR", eir_val, "bitesBed", bites_val, "dITN0", res_val, "itncov", cov_val, "Q0", Q0_val, sep = "_")]] <- list(
        split = split_name,
        init_EIR = eir_val,
        bites_Bed = bites_val,
        d_ITN0 = res_val,
        itn_cov = cov_val,
        Q0 = Q0_val,
        best_fit = best_fit,
        error = error[index]
      )
    }
    }
      }
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
    d_ITN0 = results[[name]]$d_ITN0,
    itn_cov = results[[name]]$itn_cov,
    Q0 = results[[name]]$Q0,
    endec_mu = results[[name]]$best_fit$endec_mu,
    wane_endec = results[[name]]$best_fit$wane_endec,
    error = results[[name]]$error
  )
}))

# View the results
print(results_df)

heatmap_data <- results_df %>%
  select(split, init_EIR, bites_Bed, d_ITN0, itn_cov, endec_mu, wane_endec) %>%
  pivot_longer(cols = c(endec_mu, wane_endec),
               names_to = "measure",
               values_to = "value")

#head(results_df_long)

ggplot(heatmap_data, aes(x = factor(init_EIR), y = factor(bites_Bed), fill = value))+
  geom_tile()+
  facet_grid(measure ~ split + d_ITN0 + itn_cov, scales = "free", labeller = label_both)+
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50")+
  labs(fill = "value")

require(ggpattern)

results_df <- results_df %>%
  mutate(d_ITN0 = round(d_ITN0, 2))



ggplot(results_df, aes(x = factor(Q0),
                       y = factor(split),
                       fill = endec_mu,
                       pattern = as.factor(wane_endec)))+
  geom_tile_pattern(color = "black", size = 0.5, pattern_density = 0.4,
                    pattern_fill = "grey", pattern_spacing = 0.1, pattern_colour = "white")+
  facet_grid(~bites_Bed + d_ITN0 + itn_cov + init_EIR, labeller = label_both, scales = "free")+
  scale_fill_gradient(low = "blue", high = "red", name = "endec_mu")+
  scale_pattern_manual(values = c("none", "stripe"))+ #ensure borders of legend vis) +
  theme_minimal()+
  labs(y="endec_cov")+
  theme(panel.spacing = unit(0.1, "lines"))+
  guides(
    pattern = guide_legend(
      override.aes = list(
        fill = "white",     # Change background color of the legend box here
        color = "grey",
        pattern = c("none", "stripe"),
        pattern_spacing = 0.01
      )
    )
  )

results_df_low_endec_cov <- results_df %>%
  filter(split == "low_cov")

results_df_high_endec_cov <- results_df %>%
  filter(split == "high_cov")

param_space_low_endec <- ggplot(results_df_low_endec_cov, aes(x = Q0, y = itn_cov, fill = endec_mu, size = as.factor(wane_endec),
                       shape = as.factor(init_EIR), alpha = 0.2))+
  geom_point()+
  geom_jitter(width = 0.03)+
  facet_grid(bites_Bed ~ d_ITN0, labeller = label_both)+
  guides(fill = "none",
         size = "none",
         alpha = "none")+
  scale_shape_manual(values = c(21, 24))+
  ggtitle("Low endectocide coverage")+
  labs(size = "endec_wane",
       fill = "endec_mu",
       shape = "init_EIR")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")+
  theme_bw()

param_space_high_endec <- ggplot(results_df_high_endec_cov, aes(x = Q0, y = itn_cov, fill = endec_mu, size = as.factor(wane_endec),
                                     shape = as.factor(init_EIR), alpha = 0.2))+
  geom_point()+
  geom_jitter(width = 0.03)+
  facet_grid(bites_Bed ~ d_ITN0, labeller = label_both)+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")+
  scale_shape_manual(values = c(21, 24))+
  ggtitle("High endectocide coverage")+
  labs(size = "endec_wane",
       fill = "endec_mu",
       shape = "init_EIR")+
  guides(alpha = "none",
         shape = "none")+
  theme_bw()

parameter_space_plots <- cowplot::plot_grid(param_space_low_endec, param_space_high_endec)
ggsave(parameter_space_plots, file = "analysis/exploring_interactions/malariasim-odin/param_space_plot_mu_wane.png")

#split this again by EIR

results_df_map <- results_df %>%
  mutate(ivm_cov = paste0("ivm_cov:", split),
         init_EIR = paste0("init_EIR:", init_EIR),
         bites_Bed = paste0("bites_Bed:", bites_Bed),
         d_ITN0 = paste0("d_ITN0:", d_ITN0),
         itn_cov = paste0("itn_cov:", itn_cov),
         Q0 = paste0("Q0:", Q0)) %>%
  select(ivm_cov, init_EIR, bites_Bed, d_ITN0, itn_cov, Q0, endec_mu, wane_endec)

results_df_map_long <- results_df_map %>%
  pivot_longer(cols = -c(endec_mu, wane_endec),names_to = "parameter", values_to = "values")


ggplot(results_df_map_long, aes(x = values, y = values, fill = endec_mu))+
  geom_tile()+
  scale_fill_gradient(low = "blue", high = "red", name = "endec_mu")

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
endec_mu_vec <- seq(0, 1, 0.01)
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
wanes <- unique(df_mod4$wane.x)

df_mod4_out <- df_mod4 %>%
  filter(init_EIR == 100 & endec_mu.y == 0.45 & wane.y == wanes[4])
df_mod3_out <- df_mod3 %>%
  filter(ivm_cov.x == 0.9, init_EIR == 100)

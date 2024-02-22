#fig-2-model-comparisons
require(tidyverse)
require(cowplot)
df_antag_ITN_IVM <- read.csv("analysis/exploring_interactions/model_output/antag_pyr_LLIN_IVM.csv") %>%
  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
                             TRUE ~ NA_character_))
df_antag_ITN <- read.csv("analysis/exploring_interactions/model_output/antag_pyr_LLIN.csv")%>%
  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
                             TRUE ~ NA_character_))
df_add_ITN_IVM <- read.csv("analysis/exploring_interactions/model_output/add_pyr_LLIN_IVM.csv")%>%
  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
                             TRUE ~ NA_character_))
df_add_ITN <- read.csv("analysis/exploring_interactions/model_output/add_pyr_LLIN.csv")%>%
  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
                             TRUE ~ NA_character_))

EIR_vals <- df_antag_ITN_IVM %>%
  group_by(species, itn_cov, d_ITN0) %>% #changed species, itn coverage and resistance
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

df_antag_ITN_IVM <- df_antag_ITN_IVM %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))

df_antag_ITN_IVM$init_EIR <- fillNAgaps(df_antag_ITN_IVM$init_EIR)

df_antag_ITN <- df_antag_ITN %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))

df_antag_ITN$init_EIR <- fillNAgaps(df_antag_ITN$init_EIR)

df_add_ITN_IVM <- df_add_ITN_IVM %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))
df_add_ITN_IVM$init_EIR <- fillNAgaps(df_add_ITN_IVM$init_EIR)

df_add_ITN <- df_add_ITN %>%
  mutate(init_EIR = case_when(t == 1 &EIR_tot == EIR_vals2[1] ~ "low",  #2
                              t == 1 &EIR_tot == EIR_vals2[2] ~ "medium", #25
                              t == 1 &EIR_tot == EIR_vals2[3] ~ "high",  #100
                              TRUE ~ NA_character_))
df_add_ITN$init_EIR <- fillNAgaps(df_add_ITN$init_EIR)

#check
unique(df_antag_ITN_IVM$init_EIR)
unique(df_antag_ITN$init_EIR)
unique(df_add_ITN_IVM$init_EIR)
unique(df_add_ITN$init_EIR)

#list_all <- as.list(df_antag_ITN_IVM, df_antag_ITN, df_add_ITN_IVM, df_add_ITN)
#df_all <- as.data.frame(rbind("do.call", list_all))
#
#head(df_all)

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
baseline_vals_antag <- df_antag_ITN %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  group_by(species, itn_cov, d_ITN0, init_EIR) %>%
  summarise(mean_prev_baseline = mean(slide_prev0to5),
         mean_EIR_baseline = mean(EIR_tot),
         mean_avhc_baseline = mean(avhc),
         mean_mu_baseline = mean(mu))

df_antag_eff <- df_antag_ITN_IVM %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  group_by(species, itn_cov, d_ITN0, init_EIR) %>%
  summarise(mean_prev = mean(slide_prev0to5),
            mean_EIR = mean(EIR_tot),
            mean_avhc = mean(avhc),
            mean_mu = mean(mu)) %>%
  left_join(baseline_vals_antag) %>%
  mutate(eff_prev_antag = mean_prev_baseline - mean_prev,
         eff_eir_antag = mean_EIR_baseline - mean_EIR)
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
#generally following a straight line :)

lm_prev <- lm(eff_prev_add ~ eff_prev_antag, data = df_eff)
summary(lm_prev) #R2 = 0.9987


pals_res <-rev(c('#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026'))

plot_prev <- ggplot(df_eff, aes(x = eff_prev_antag, eff_prev_add))+
  geom_point()+
  ylab("Absolute reduction in prevalence (additive model)")+
  xlab("Absolute reduction in prevalence (antagonistic model)")+
  theme_minimal()+
  stat_smooth(method = "lm", col = "red", linetype = "dashed")+
  annotate("text", label = "italic(R)^2 == 0.9987", parse = TRUE, x =0.04, y = 0.06, col = "red")

#helps to see that some clustering  but generally following a straight line
#change this so appears in order: low, medium, high
plot_prev_facet <- ggplot(df_eff, aes(x = eff_prev_antag, eff_prev_add, shape = as.factor(species), col = as.factor(d_ITN0)), stroke = 1)+
  geom_point()+
  facet_grid(rows = vars(itn_cov), cols = vars(init_EIR))+
  ylab("Absolute reduction in prevalence (additive model)")+
  xlab("Absolute reduction in prevalence (antagonistic model)") +
  theme_minimal()+
  labs(shape = "Anopheles species")+
  scale_colour_manual(name = "Resistance (% survival in bioassay)",
                      labels = c("90%", "70%", "50%", "10%", "0%, no resistance"),
                      values = pals_res)


#for eir

lm_eir <- lm(eff_eir_add ~ eff_eir_antag, data = df_eff)
summary(lm_eir) #R2 =  0.9995
plot_EIR <- ggplot(df_eff, aes(x = eff_eir_antag, eff_eir_add))+
  geom_point()+
  ylab("Absolute reduction in EIR (additive model)")+
  xlab("Absolute reduction in EIR (antagonistic model)")+
  stat_smooth(method = "lm", col = "red", linetype = "dashed")+
  annotate("text", label = "italic(R)^2 == 0.9995", parse = TRUE, x =50, y = 75, col = "red")+
  theme_minimal()


#helps to see that some clustering  but generally following a straight line
#change this so appears in order: low, medium, high
plot_EIR_facet <- ggplot(df_eff, aes(x = eff_eir_antag, eff_eir_add, shape = as.factor(species), col = as.factor(d_ITN0)))+
  geom_point()+
  facet_grid(rows = vars(itn_cov), cols = vars(init_EIR))+
  ylab("Absolute reduction in EIR (additive model)")+
  xlab("Absolute reduction in EIR (antagonistic model)") +
  theme_minimal()+
  labs(shape = "Anopheles species") #+
  scale_colour_manual(name = "Resistance (% survival in bioassay)",
                      labels = c("90%", "70%", "50%", "10%", "0%, no resistance"),
                      values = pals_res)
#
model_comparisons_plot <- plot_grid(plot_prev, plot_EIR, ncol = 2)

model_comparisons_plot_prev <- plot_grid(plot_prev, plot_prev_facet)

model_comparisons_plot_eir <- plot_grid(plot_EIR, plot_EIR_facet)

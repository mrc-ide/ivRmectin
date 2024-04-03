#fig-2-model-comparisons
require(tidyverse)
require(cowplot)
df_antag_ITN_IVM <- read_rds("analysis/exploring_interactions/model_output/antag_pyr_LLIN_IVM_nets_age.rds") %>%
  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
                             TRUE ~ NA_character_))
df_antag_ITN <- read_rds("analysis/exploring_interactions/model_output/antag_pyr_LLIN_nets_age.rds")%>%
  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
                             TRUE ~ NA_character_))
df_add_ITN_IVM <- read_rds("analysis/exploring_interactions/model_output/add_pyr_LLIN_IVM.rds")%>%
  mutate(species = case_when(bites_Bed == 0.85 & Q0 == 0.92 ~ "gambiae",
                             bites_Bed == 0.8 & Q0 == 0.71 ~ "arabiensis",
                             bites_Bed == 0.78 & Q0 == 0.94 ~ "funestus",
                             bites_Bed == 0.52 & Q0 == 0.21 ~ "stephensi",
                             TRUE ~ NA_character_))
df_add_ITN <- read_rds("analysis/exploring_interactions/model_output/add_pyr_LLIN.rds")%>%
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

##filter to medium EIR, 10% resistance and 80% LLIN coverage

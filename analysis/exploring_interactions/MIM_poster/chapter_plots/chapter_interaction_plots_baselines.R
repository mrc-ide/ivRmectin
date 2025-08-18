#we have counterfactual of ivermectin only or no intervention

#ivm only
antag_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/Q0/antag_IVM.rds")

#no intervention
antag_no_int <- readRDS("analysis/exploring_interactions/MIM_poster/Q0/antag_no_int.rds")

#bites bed####
bb_antag <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/antag.rds")

bb_antag_LLIN_IVM <- bb_antag %>%
  filter(model == "antag_LLIN_IVM")
bb_antag_LLIN <- bb_antag %>%
  filter(model == "antag_LLIN")

#d_ITN0
dITN0_antag <- readRDS("analysis/exploring_interactions/MIM_poster/itn_res/antag.rds")

dITN0_antag_LLIN_IVM <- dITN0_antag %>%
  filter(model == "antag_LLIN_IVM")
dITN0_antag_LLIN <- dITN0_antag %>%
  filter(model == "antag_LLIN")

#itn_cov
itn_cov_antag <- readRDS("analysis/exploring_interactions/MIM_poster/itn_cov/antag.rds")

itn_cov_antag_LLIN_IVM <- itn_cov_antag %>%
  filter(model == "antag_LLIN_IVM")
itn_cov_antag_LLIN <- itn_cov_antag %>%
  filter(model == "antag_LLIN")

#Q0
Q0_antag <- readRDS("analysis/exploring_interactions/MIM_poster/Q0/antag.rds")

Q0_antag_LLIN_IVM <- Q0_antag %>%
  filter(model == "antag_LLIN_IVM")
Q0_antag_LLIN <- Q0_antag %>%
  filter(model == "antag_LLIN")



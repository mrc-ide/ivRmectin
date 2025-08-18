require(tidyverse)

df_bb_Q0_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN_IVM.rds")
df_bb_cov_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN_IVM.rds")
df_bb_res_ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN_IVM.rds")

df_bb_Q0_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_Q0_ITN.rds")
df_bb_cov_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN.rds")
df_bb_res_ITN <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_res_ITN.rds")

net_seq <- seq(100, 3650, by = 3*365)
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution

#heatmap bites_Bed and Q0

output_bb_Q0_ITN_IVM <- df_bb_Q0_ITN_IVM %>%
  select(t, bites_Bed, Q0, clin_inc0to5) %>%
  rename(clin_inc_bb_Q0_ITN_IVM = clin_inc0to5)


output_bb_Q0_ITN <- df_bb_Q0_ITN %>%
  select(t, bites_Bed, Q0, clin_inc0to5) %>%
  rename(clin_inc_ITN = clin_inc0to5)

output_bb_Q0 <- left_join(output_bb_Q0_ITN, output_bb_Q0_ITN_IVM)

output_bb_Q0 <- output_bb_Q0 %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(Q0, bites_Bed) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_bb_Q0_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100)

ggplot(output_bb_Q0, aes(x = bites_Bed, y = Q0, fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()

#heatmap bites_Bed and cov

output_bb_cov_ITN_IVM <- df_bb_cov_ITN_IVM %>%
  select(t, bites_Bed, itn_cov, clin_inc0to5) %>%
  rename(clin_inc_bb_cov_ITN_IVM = clin_inc0to5)


output_bb_cov_ITN <- df_bb_cov_ITN %>%
  select(t, bites_Bed, itn_cov, clin_inc0to5) %>%
  rename(clin_inc_ITN = clin_inc0to5)

output_bb_cov <- left_join(output_bb_cov_ITN, output_bb_cov_ITN_IVM)

output_bb_cov <- output_bb_cov %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(bites_Bed, itn_cov) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_bb_cov_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100)

ggplot(output_bb_cov, aes(x = bites_Bed, y = itn_cov, fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()

#heatmap bites_bed and res
output_bb_res_ITN_IVM <- df_bb_res_ITN_IVM %>%
  select(t, bites_Bed, d_ITN0, r_ITN0,  clin_inc0to5) %>%
  rename(clin_inc_bb_res_ITN_IVM = clin_inc0to5)


output_bb_res_ITN <- df_bb_res_ITN %>%
  select(t, bites_Bed, d_ITN0, r_ITN0, clin_inc0to5) %>%
  rename(clin_inc_ITN = clin_inc0to5)

output_bb_res <- left_join(output_bb_res_ITN, output_bb_res_ITN_IVM)

output_bb_res <- output_bb_res %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(bites_Bed, d_ITN0, r_ITN0) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_bb_res_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100)

ggplot(output_bb_res, aes(x = bites_Bed, y = d_ITN0, fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()
##with high bites_Bed, there is little variability in impact due to d_ITN0. Check inputs

ggplot(output_bb_res, aes(x = bites_Bed, y = r_ITN0, fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw()

require(tidyverse)
#heatmap for efficacy across different Q0 and IVM scenarios
ITN_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/Q0_ITN_IVM_cov.rds")
ITN_only <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/Q0_ITN_only.rds")

output_ITN_IVM <- ITN_IVM %>%
  select(t, ivm_cov, Q0, clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  rename(clin_inc_Q0_ITN_IVM = clin_inc0to5)


output_ITN <- ITN_only %>%
  select(t, Q0, clin_inc0to5) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>%
  rename(clin_inc_ITN = clin_inc0to5)

output_Q0 <- left_join(output_ITN, output_ITN_IVM)
net_seq <- seq(100, 3650, by = 3*365)
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution


output_Q0 <- output_Q0 %>%
  filter(between(t, IVM_begin1, IVM_begin1+180)) %>%
  group_by(Q0, ivm_cov) %>%
  summarise(tot_cases_ITN_IVM = sum(clin_inc_Q0_ITN_IVM),
            tot_cases_ITN = sum(clin_inc_ITN)) %>%
  mutate(abs_diff_IVM = tot_cases_ITN - tot_cases_ITN_IVM,
         rel_diff_IVM = ((tot_cases_ITN - tot_cases_ITN_IVM)/tot_cases_ITN)*100)

ivm_cov_el <- unique(output_Q0$ivm_cov)
ivm_cov_in <- c(0.3, 0.5, 0.7, 0.9)

output_Q0 <- output_Q0 %>%
  mutate(ivm_cov = case_when(ivm_cov == ivm_cov_el[1] ~ ivm_cov_in[1],
                             ivm_cov == ivm_cov_el[2] ~ ivm_cov_in[2],
                             ivm_cov == ivm_cov_el[3] ~ ivm_cov_in[3],
                             TRUE ~ ivm_cov_in))

common_limits_rel <- c(10, 45)
common_limits_abs <- c(50,600)


heatmap_Q0_ivm_rel <- ggplot(output_Q0, aes(x = as.factor(Q0), y = as.factor(ivm_cov*100), fill = rel_diff_IVM))+
  geom_tile()+
  theme_bw(base_size = 16)+
  scale_fill_viridis_c(limits = common_limits_rel, name = "Cases averted (%) in \n under 5-year-olds due to endectocide")+
  xlab("Human Blood Index")+
  ylab("Ivermectin coverage (%)")+
  #guides(fill = "none")+
  #geom_text(aes(label = species), parse = TRUE, col = "white", size = 5)+
  geom_text(aes(label = paste0(round(rel_diff_IVM, 1), "%")),
            col = "white", size = 5)+
  theme(legend.position = "bottom")
ggsave(heatmap_Q0_ivm_rel, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/heatmap_Q0_ivm_cov_rel.pdf")


heatmap_Q0_ivm_abs <- ggplot(output_Q0, aes(x = as.factor(Q0), y = as.factor(ivm_cov*100), fill = abs_diff_IVM))+
  geom_tile()+
  theme_bw(base_size = 16)+
  scale_fill_fermenter(limits = common_limits_abs, name = "Absolute reduction in cases in \n under 5-year-olds due to endectocide")+
  xlab("Human Blood Index")+
  ylab("Ivermectin coverage (%)")+
  #guides(fill = "none")+
  #geom_text(aes(label = species), parse = TRUE, col = "white", size = 5)+
  geom_text(aes(label = paste0(round(abs_diff_IVM))),
            col = "black", size = 5)+
  theme(legend.position = "bottom",
        legend.direction = "vertical")

heatmap_Q0_ivm <- cowplot::plot_grid(heatmap_Q0_ivm_rel, heatmap_Q0_ivm_abs,
                                     labels = c("A", "B"))

ggsave(heatmap_Q0_ivm, file = "analysis/exploring_interactions/MIM_poster/chapter_plots/heatmap_Q0_ivm_cov.pdf")

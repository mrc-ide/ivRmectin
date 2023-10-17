#efficacy of LLIN and IVM combo, LLIN only compared to no interventions
#under different coverage and resistance scenarios

require(tidyverse)

#choose palette based on interventions

pals <- c( "#DDCC77", "#117733", "#332288", "#AA4499",
          "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


#breaks_plot <- seq(2250, 3650, by = 30)
breaks <- seq(2250, 3650, by = 30*6)
time_period <- 365*10
IVM_begin <- (365*8)+180
mda_int <- 30
IVM_start <- c(IVM_begin, IVM_begin+mda_int, IVM_begin+mda_int+mda_int)

eff_len <- 23
ivm_on <- IVM_start[1] #3100
ivm_off <- IVM_start[3]+eff_len #3183

net_seq <- seq(100, 3650, by = 3*365)


spikes <-  read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE)
res_read <- read.csv("analysis/esa-analysis/output/pyr_out_antag_IVM.csv", header = TRUE)
res <- unique(res_read$d_ITN0)
rm1 <- spikes %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2]) %>% #10% res
  filter(between(t, 1190, 1199)) #remove t = 1195, t = 1196, t = 1197

spikes %>%
  filter(itn_cov == 0.8 & d_ITN0 ==  res[2]) %>%
  filter(between(t, 3375, 3395)) #remove t = 2290, t = 2291, t = 2292

rm <- c(2290, 2291, 2292, 3385, 3386, 3387)


#read in no intervention file

df_no_int <- read.csv("analysis/esa-analysis/output/df_no_int_model.csv", header = TRUE) %>%
  mutate(model = "base model", int = "no int")

ggplot(df_no_int, aes(x = t, y = EIR_tot))+
  geom_line()+
  theme_minimal()

df_no_int$EIR_tot[1]

df_no_int_summary <- df_no_int %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5))

#mean EIR:  75.48611
#mean prev:  0.636088

df_no_int_summary$mean_EIR





#read in the other files
antag_mod_LLIN<- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov, EIR_tot, slide_prev0to5) %>%
  mutate(model = "antagonistic", int = "LLIN")

antag_mod_LLIN %>%
  filter(itn_cov == 0.8) %>%
  ggplot(aes(x = t, y = EIR_tot, colour = as.factor(d_ITN0)))+
  geom_line()+
  geom_line(data = df_no_int, aes(x = t, y = EIR_tot))+
  annotate(geom = "rect", xmin = ivm_on, xmax = ivm_off, ymin = 0, ymax = Inf, fill = "pink", alpha =0.2)

antag_mod_LLIN %>%
  group_by(d_ITN0, itn_cov) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(mean_EIR = mean(EIR_tot),
            mean_prev = mean(slide_prev0to5))

add_mod_LLIN<- read.csv("analysis/esa-analysis/output/pyr_out_add.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov,EIR_tot, slide_prev0to5) %>%
  mutate(model = "additive", int = "LLIN")

antag_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_antag_IVM.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov,EIR_tot, slide_prev0to5) %>%
  mutate(model = "antagonistic", int = "LLIN & IVM")

add_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_add_IVM.csv", header = TRUE) %>%
  select(t, mv, d_ITN0, itn_cov,EIR_tot, slide_prev0to5) %>%
  mutate(model = "additive", int = "LLIN & IVM")

files <- list(antag_mod_LLIN, add_mod_LLIN, antag_mod_LLIN_IVM, add_mod_LLIN_IVM)

ggplot(antag_mod_LLIN, aes(x = t, y = EIR_tot, colour = as.factor(d_ITN0)))+
  geom_line()

res_vec <- unique(antag_mod_LLIN_IVM$d_ITN0)

summary_function <- function(x) {
  res_vec <- unique(antag_mod_LLIN_IVM$d_ITN0)
  x %>%
    group_by(d_ITN0, itn_cov) %>%
    filter(between(t, ivm_on, ivm_off)) %>%
    summarise(mean_EIR = mean(EIR_tot),
              mean_prev = mean(slide_prev0to5),
              model = model,
              int = int) %>%
    mutate(res = case_when(d_ITN0 == res_vec[1] ~ "No resistance",
                         d_ITN0 == res_vec[2] ~ "10",
                         d_ITN0 == res_vec[3] ~ "50",
                         d_ITN0 == res_vec[4] ~ "70",
                         d_ITN0 == res_vec[5] ~ "90",
                         TRUE ~ NA_character_)) %>%
    filter(res != "No resistance")
}

files_summary <- files %>%
  lapply(summary_function)

itn_cov_vector <- seq(0, 0.8, 0.2)

files_summary_df  <- do.call(rbind,
                               sapply(1:(length(files_summary)), function(x){
                                 as.data.frame(files_summary[[x]]) %>%
                                   mutate(ref = x)
                               }, simplify = F))



efficacy_summary <- files_summary_df %>%
  mutate(abs_diff_EIR = df_no_int_summary$mean_EIR- mean_EIR,
         rel_diff_EIR = ((df_no_int_summary$mean_EIR - mean_EIR)/df_no_int_summary$mean_EIR)*100,
         abs_diff_prev =  df_no_int_summary$mean_prev- mean_prev,
         rel_diff_prev = ((df_no_int_summary$mean_prev - mean_prev)/df_no_int_summary$mean_prev)*100) %>%
  filter(itn_cov %in% c(0.4, 0.6, 0.8))

efficacy_summary %>%
  filter(res == "90" & itn_cov == 0.8)



efficacy_summary$plot <- factor(efficacy_summary$model, levels = c("antagonistic", "additive"),
                                labels = c("Antagonistic model", "Additive model"))

check <- efficacy_summary %>%
  filter(itn_cov == 0.8) %>%
  group_by(plot, res)

ggplot(efficacy_summary, aes(x = itn_cov, y = rel_diff_EIR, fill = as.factor(res)))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_grid(cols = vars(model), rows = vars(int))+
  scale_fill_manual(values = pals[4:1], name = "Resistance level")



#this generally behaves strangely with coverage, going to just run at 80% coverage across a range of resistance profiles


efficacy_plot_EIR <- efficacy_summary %>%
  filter(itn_cov == 0.8) %>%
  ggplot(aes(x = res, y = rel_diff_EIR, fill = as.factor(plot)))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(values = c(pals[5], pals[7]), name = "Model")+
  facet_grid(~int)+
  theme_minimal()+
  ylab("Relative difference (%) in EIR \n due to intervention")+
  xlab("Resistance (%)")+
  ylim(-20, 100)+
  theme(legend.position = "none")

efficacy_plot_prev<- efficacy_summary %>%
  filter(itn_cov == 0.8) %>%
  ggplot(aes(x = res, y = rel_diff_prev, fill = as.factor(plot)))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(values = c(pals[5], pals[7]), name = "Model")+
  facet_grid(~int)+
  theme_minimal()+
  ylab("Relative difference (%) in prevalence \n due to intervention")+
  xlab("Resistance (%)")+
  ylim(0, 80)+
  theme(legend.position = c(0.1, 0.8))


#efficacy_plot_EIR <- efficacy_summary %>%
#  filter(itn_cov == 0.8) %>%
#  ggplot(aes(x = res, y = rel_diff_EIR, fill = as.factor(plot)))+
#  geom_bar(stat = "identity", position = position_dodge())+
#  scale_fill_manual(values = c(pals[7], pals[9]), name = "model")+
#  facet_grid(~res)+
#  theme_minimal()+
#  ylab("Relative difference in EIR (%) \n compared to no interventions")+
#  xlab("Intervention")+
#  theme(legend.position = c(0.3, 0.1), legend.direction = "horizontal", text =
#        element_text(size = 12))

#efficacy_plot_prev <- efficacy_summary %>%
#  filter(itn_cov == 0.8) %>%
#  ggplot(aes(x = int, y = rel_diff_prev, fill = as.factor(res)))+
#  geom_bar(stat = "identity", position = position_dodge())+
#  scale_fill_manual(values = pals[4:1], name = "Resistance level")+
#  facet_grid(~plot)+
#  theme_minimal()+
#  ylab("Relative difference in prevalence (under 5s) (%) \n compared to no interventions")+
#  xlab("Intervention")+
#  theme(legend.position = "none", text =
#          element_text(size = 12))

require(cowplot)
efficacy_plots <- plot_grid(efficacy_plot_EIR, efficacy_plot_prev, labels = c("A)", "B)"), nrow = 2,
                            label_x = 0.1)
#ggsave(efficacy_plots, file = "plots/esa-figures/fig_3_efficacy_plot.svg")

compare_EIR_models <- efficacy_summary %>%
  filter(itn_cov == 0.8) %>%
  filter(int == "LLIN & IVM") %>%
  group_by(res,model) %>%
  summarise(rel_diff_EIR = rel_diff_EIR) %>%
  distinct() %>%
  spread(model, rel_diff_EIR) %>%
  mutate(abs_diff_EIR_model = additive-antagonistic)

compare_prev_models <- efficacy_summary %>%
  filter(itn_cov == 0.8) %>%
  filter(int == "LLIN & IVM") %>%
  group_by(res,model) %>%
  summarise(rel_diff_prev = rel_diff_prev) %>%
  distinct() %>%
  spread(model, rel_diff_prev) %>%
  mutate(abs_diff_prev_model = additive-antagonistic)



#show tom
ggplot(df_no_int, aes(x = t, y = EIR_tot))+
  geom_line()+
  theme_minimal()

antag_mod_LLIN <- read.csv("analysis/esa-analysis/output/pyr_out_antag.csv", header = TRUE)
add_mod_LLIN <- read.csv("analysis/esa-analysis/output/pyr_out_add.csv", header = TRUE)
antag_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_antag_IVM.csv", header = TRUE)
add_mod_LLIN_IVM <- read.csv("analysis/esa-analysis/output/pyr_out_add_IVM.csv", header = TRUE)

#####################################

antag_LLIN <- antag_mod_LLIN %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
            mean_prev = mean(slide_prev0to5)) %>%
  mutate(int = "LLIN")

antag_both <- antag_mod_LLIN_IVM %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
            mean_prev = mean(slide_prev0to5)) %>%
  mutate(int = "LLIN & IVM")

summary_antag <- rbind(antag_LLIN, antag_both) %>%
  select(-model)
summary_wide_antag <- spread(summary_antag, int, mean_prev) %>%
  mutate(diff = LLIN - `LLIN & IVM`,
         rel_diff = diff/LLIN)


add_LLIN <- add_mod_LLIN %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
    mean_prev = mean(slide_prev0to5)) %>%
  mutate(int = "LLIN")

add_both <- add_mod_LLIN_IVM %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
    mean_prev = mean(slide_prev0to5)) %>%
  mutate(int = "LLIN & IVM")

summary_add <- rbind(add_LLIN, add_both) %>%
  select(-model)
summary_wide_add <- spread(summary_add, int, mean_prev) %>%
  mutate(diff = LLIN - `LLIN & IVM`,
         rel_diff = diff/LLIN)

#summary wide antag and summary wide add show the additional reduction in prevalence that ivermectin gives on top of LLINs, across resistance profiles
#we see that across both models, in terms of resource allocation, it would be more cost-effective to distribute ivermectin to an area of higher insecticide resistance as it mops up more residual transmission

######

antag_LLIN_EIR <- antag_mod_LLIN %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
    mean_EIR = mean(EIR_tot)) %>%
  mutate(int = "LLIN")

antag_both_EIR <- antag_mod_LLIN_IVM %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
    mean_EIR = mean(EIR_tot)) %>%
  mutate(int = "LLIN & IVM")

summary_antag_EIR<- rbind(antag_LLIN_EIR, antag_both_EIR) %>%
  select(-model)
summary_wide_antag_EIR <- spread(summary_antag_EIR, int, mean_EIR) %>%
  mutate(diff = LLIN - `LLIN & IVM`,
         rel_diff = diff/LLIN)


add_LLIN_EIR <- add_mod_LLIN %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
    mean_EIR = mean(EIR_tot)) %>%
  mutate(int = "LLIN")

add_both_EIR <- add_mod_LLIN_IVM %>%
  filter(itn_cov == 0.8) %>%
  group_by(d_ITN0, model) %>%
  filter(between(t, ivm_on, ivm_off)) %>%
  summarise(
    mean_EIR = mean(EIR_tot)) %>%
  mutate(int = "LLIN & IVM")

summary_add_EIR <- rbind(add_LLIN_EIR, add_both_EIR) %>%
  select(-model)
summary_wide_add_EIR <- spread(summary_add_EIR, int, mean_EIR) %>%
  mutate(diff = LLIN - `LLIN & IVM`,
         rel_diff = diff/LLIN)


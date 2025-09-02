require(tidyverse)


#initial setup
itn_on <- 100 #introduce nets 100 days into simulation

net_seq <- seq(100, 3650, by = 3*365)
mda_int <- 30
#ivm on when nets are 6 months old
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)


#checking dynamics and trends etc in bites_Bed and influence on IVM impact
df_bb_cov_ITN_IVM2 <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN_IVM.rds")
df_bb_cov_ITN2 <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/bb_cov_ITN.rds")

df_bb <- rbind(df_bb_cov_ITN_IVM2, df_bb_cov_ITN2) %>%
  filter(itn_cov == 0.8)

ggplot(df_bb, aes(x = t, y = EIRout*1000, lty = as.factor(int), col = as.factor(bites_Bed)))+
  geom_line()

ggplot(df_bb, aes(x = t, y = clin_inc0to5*1000, lty = as.factor(int), col = as.factor(bites_Bed)))+
  geom_line()

bb_IVM_cov80_dynamics <-ggplot(df_bb, aes(x = t, y = slide_prev0to5*100, lty = as.factor(int), col = as.factor(bites_Bed)))+
  geom_line()+
  theme_bw()+
  ylim(0, 70)

df_bb_all <- rbind(df_bb_cov_ITN_IVM2, df_bb_cov_ITN2)

ggplot(df_bb_all, aes(x = t, y = EIRout*1000, lty = as.factor(int), col = as.factor(bites_Bed)))+
  geom_line()+
  facet_wrap(vars(itn_cov), labeller = label_both)+
  geom_vline(xintercept = IVM_begin1)+
  geom_vline(xintercept = IVM_begin1+180)

ggplot(df_bb_all, aes(x = t, y = clin_inc0to5*1000, lty = as.factor(int), col = as.factor(bites_Bed)))+
  geom_line()+
  facet_wrap(vars(itn_cov), labeller = label_both)+
  geom_vline(xintercept = IVM_begin1)+
  geom_vline(xintercept = IVM_begin1+180)

prev_dynamics_plot <- ggplot(df_bb_all, aes(x = t, y = slide_prev0to5*100, lty = as.factor(int), col = as.factor(bites_Bed)))+
  geom_line()+
  facet_wrap(vars(itn_cov), labeller = label_both)+
  geom_vline(xintercept = IVM_begin1)+
  geom_vline(xintercept = IVM_begin1+180)+
  ylim(0,100)+
  theme_bw()

#the absolute reduction for all of these is bigger if there is more residual transmission

out_bb_cov_ITN_IVM_prev<- df_bb_cov_ITN_IVM2 %>%
  select(t, bites_Bed, itn_cov,slide_prev0to5) %>%
  filter(t == IVM_start1[3]+30) %>%
  group_by(itn_cov, bites_Bed) %>%
  summarise(prev_bb_ITN_IVM = slide_prev0to5)

out_bb_cov_ITN_prev<- df_bb_cov_ITN2 %>%
  select(t, bites_Bed, itn_cov,slide_prev0to5) %>%
  filter(t == IVM_start1[3]+30) %>%
  group_by(itn_cov, bites_Bed) %>%
  summarise(prev_bb_ITN = slide_prev0to5)

out_bb_cov <- left_join(out_bb_cov_ITN_prev, out_bb_cov_ITN_IVM_prev)

out_bb_cov_comp <- out_bb_cov %>%
  mutate(abs_diff_IVM = prev_bb_ITN - prev_bb_ITN_IVM,
         rel_diff_IVM = ((prev_bb_ITN - prev_bb_ITN_IVM)/prev_bb_ITN)*100)

#for high ITN cov, with lowest bb, the abs reduction due to ivm is 0.127, compared to 0.0363 at high bb level
#this is a 70% difference in predicted efficacy, influenced by bites_Bed

#MATMAL was powered to detect a 50% difference (10% in control, 5% in int)
#but for example in a high transmission setting with high net cov and bites_Bed and low res, only going to see an 20% reduction due to ivermectin

#for cases averted
out_bb_cov_ITN_IVM_inc<- df_bb_cov_ITN_IVM2 %>%
  select(t, bites_Bed, itn_cov,clin_inc0to5, EIRout) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>% #make it per 1000 persons
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  group_by(itn_cov, bites_Bed) %>%
  summarise(inc_bb_ITN_IVM = sum(clin_inc0to5),
            EIR_bb_ITN_IVM = sum(EIRout))

out_bb_cov_ITN_inc<- df_bb_cov_ITN2 %>%
  select(t, bites_Bed, itn_cov,clin_inc0to5, EIRout) %>%
  mutate(clin_inc0to5 = clin_inc0to5*1000) %>% #make it per 1000 persons
  filter(between(t, IVM_start1[1],IVM_start1[1]+180)) %>%
  group_by(itn_cov, bites_Bed) %>%
  summarise(inc_bb_ITN = sum(clin_inc0to5),
            EIR_bb_ITN = sum(EIRout))

out_bb_cov_inc <- left_join(out_bb_cov_ITN_inc, out_bb_cov_ITN_IVM_inc)

out_bb_cov_comp_inc <- out_bb_cov_inc %>%
  mutate(inc_abs_diff_IVM = inc_bb_ITN - inc_bb_ITN_IVM,
         inc_rel_diff_IVM = ((inc_bb_ITN - inc_bb_ITN_IVM)/inc_bb_ITN)*100,
         EIR_abs_diff_IVM = EIR_bb_ITN - EIR_bb_ITN_IVM,
         EIR_rel_diff_IVM = ((EIR_bb_ITN - EIR_bb_ITN_IVM)/EIR_bb_ITN)*100)

heatmap_bb_itn_cov_rel <- ggplot(out_bb_cov_comp_inc,
       aes(x = as.factor(bites_Bed),
           y = as.factor(itn_cov*100),
           fill = inc_rel_diff_IVM)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cases averted (%) in \n under 5-year-olds due to endectocide") +
  geom_text(aes(label = paste0(round(inc_rel_diff_IVM, 1), "%")),
            color = "white") +
 labs(y = "ITN coverage (%)", x = "Proportion of bites when people are in bed")+
  theme_bw()

heatmap_bb_itn_cov_abs <- ggplot(out_bb_cov_comp_inc,
                                 aes(x = as.factor(bites_Bed),
                                     y = as.factor(itn_cov*100),
                                     fill = inc_abs_diff_IVM)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cases averted (absolute) in \n under 5-year-olds per 1000 persons due to endectocide") +
  geom_text(aes(label = round(inc_abs_diff_IVM, 3)),
            color = "white")+
  labs(y = "ITN coverage (%)", x = "Proportion of bites when people are in bed")

heatmap_bb_itn_cov_rel_EIR <- ggplot(out_bb_cov_comp_inc,
                                 aes(x = as.factor(bites_Bed),
                                     y = as.factor(itn_cov*100),
                                     fill = EIR_rel_diff_IVM)) +
  geom_tile() +
  scale_fill_viridis_c(name = "EIR averted (%) in \n under 5-year-olds due to endectocide") +
  geom_text(aes(label = paste0(round(EIR_rel_diff_IVM, 1), "%")),
            color = "white")+
  labs(y = "ITN coverage (%)", x = "Proportion of bites when people are in bed")+
  theme_bw()

###########SHOW THEM THIS##################
heatmap_bb_itn_cov_abs_EIR <- ggplot(out_bb_cov_comp_inc,
                                     aes(x = as.factor(bites_Bed),
                                         y = as.factor(itn_cov*100),
                                         fill = EIR_abs_diff_IVM)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", name = "EIR averted (absolute) in \n under 5-year-olds due to endectocide") +
  geom_text(aes(label = round(EIR_abs_diff_IVM, 1)),
            color = "white")+
  labs(y = "ITN coverage (%)", x = "Proportion of bites when people are in bed")+
  theme_bw()
############################################

#re-run this without ITNs (so baseline no interventions vs ivermectin only)
#what is the predicted impact?

Q0_baseline <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/Q0_no_int.rds")
Q0_IVM <- readRDS("analysis/exploring_interactions/MIM_poster/bites_Bed/Q0_IVM.rds") #IVM coverage is 0.7

Q0_df <- rbind(Q0_baseline, Q0_IVM)

ggplot(Q0_df, aes(x = t, y = EIRout*1000, lty = as.factor(int), col = as.factor(Q0)))+
  geom_line()+
  ylim(0, 30)

#for cases averted
out_Q0_IVM <- Q0_IVM %>%
  select(t, Q0, clin_inc0to5, EIRout) %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  group_by(Q0) %>%
  summarise(inc_Q0_IVM = sum(clin_inc0to5),
            EIR_Q0_IVM = sum(EIRout))

out_Q0_baseline <- Q0_baseline %>%
  select(t, Q0, clin_inc0to5, EIRout) %>%
  filter(between(t, IVM_start1[1], IVM_start1[1]+180)) %>%
  group_by(Q0) %>%
  summarise(inc_Q0 = sum(clin_inc0to5),
            EIR_Q0 = sum(EIRout))

out_Q0 <- left_join(out_Q0_IVM, out_Q0_baseline)

out_Q0_comp <- out_Q0 %>%
  mutate(inc_abs_diff_IVM = inc_Q0 - inc_Q0_IVM,
         inc_rel_diff_IVM = ((inc_Q0 - inc_Q0_IVM)/inc_Q0)*100,
         EIR_abs_diff_IVM = EIR_Q0 - EIR_Q0_IVM,
         EIR_rel_diff_IVM = ((EIR_Q0 - EIR_Q0_IVM)/EIR_Q0)*100)

#how much would a 48.5% vs 40.4% cases averted estimate influence power calculations?

ggplot(out_Q0_comp, aes(x = as.factor(Q0), y = inc_rel_diff_IVM ))+
  geom_bar(stat = "identity")

ggplot(out_Q0_comp, aes(x = as.factor(Q0), y = EIR_rel_diff_IVM ))+
  geom_bar(stat = "identity")

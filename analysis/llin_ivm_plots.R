require(tidyverse)
require(cowplot)
#av_mosq and avhc under different LLIN coverages
out.df <- read.csv("data/out_df.csv", header = TRUE) #made in llin_ivm.R

#PLOT 1
#plot of av_mosq_sum (across the intervention categories changing with ITN cov). Difficult to interpret because summing
ggplot(out.df, aes(x = t, y = av_mosq_sum, col = as.factor(itn_cov)))+
  geom_line(linewidth = 1)

#PLOT 2
#plot of avhc over time - this is the normalised value and what clearly affects the ivm uptake
avhc_plot <- ggplot(out.df, aes(x = t, y = avhc, col = as.factor(itn_cov)))+
  geom_line(linewidth = 1)+
  ylim(0, 1)


mu_plot <- ggplot(out.df, aes(x = t, y = mu, col = as.factor(itn_cov)))+
  geom_line(linewidth = 1)+
  ylim(0, 1)


plot_grid(avhc_plot, mu_plot, ncol = 2)

high_endo <- out.df2 %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc)) %>%
  ggplot()+
  geom_point(aes(x = itn_cov, y = mean_avhc))

#looking at how uptake of ivermectin by mosquitoes changes with LLIN cov and IVM cov, high endophagy
out.df2 <- read.csv("data/out_df_2.csv", header = TRUE)

high_endo <- out.df2 %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc)) %>%
  ggplot()+
  geom_point(aes(x = itn_cov, y = mean_avhc))+
  ylim(0, 0.5)

checking <- out.df2 %>%
  filter(Extot <0 | Ixtot <0|Sxtot < 0)
unique(checking$ivm_cov)
unique(checking$itn_cov)

range(checking$Sxtot)
range(checking$Extot)
range(checking$Ixtot)

unique(out.df2$itn_cov)
unique(out.df2$ivm_cov)

out.df2 <- out.df2 %>%
  mutate(is_neg  = case_when(if_any(contains("xtot"), ~. <0) ~ "yes",
                             TRUE ~ "no"),
         prop_Sxtot = Sxtot/mvxtot,
         prop_Extot = Extot/mvxtot,
         prop_Ixtot = Ixtot/mvxtot)


show <- out.df2 %>%
  filter(itn_cov == 1 & is_neg == "yes")

head(show)

#low bites_Bed and bites_Indoors file, looping through different IVM and LLIN coverages
out.df3 <- read.csv("data/out_df_3.csv", header = TRUE)

low_endo <- out.df3 %>%
  group_by(itn_cov) %>%
  summarise(mean_avhc = mean(avhc)) %>%
  ggplot()+
  geom_point(aes(x = itn_cov, y = mean_avhc))+
  ylim(0, 0.5)


avhc_endophagy <- cowplot::plot_grid(high_endo, low_endo)


out.df3 <- out.df3 %>%
  mutate(is_neg  = case_when(if_any(contains("xtot"), ~. <0) ~ "yes",
                             TRUE ~ "no"),
         prop_Sxtot = Sxtot/mvxtot,
         prop_Extot = Extot/mvxtot,
         prop_Ixtot = Ixtot/mvxtot)

#can put this in a plotting script
Sx_plot <- ggplot(out.df2, aes(x = t, y = Sxtot, col = as.factor(is_neg)))+
  geom_point()+
  geom_line()+
  facet_grid(itn_cov ~ ivm_cov)

Ex_plot <- ggplot(out.df2, aes(x = t, y = Extot, col = as.factor(is_neg)))+
  geom_point()+
  geom_line()+
  facet_grid(itn_cov ~ ivm_cov)

Ix_plot <- ggplot(out.df2, aes(x = t, y = Ixtot, col = as.factor(is_neg)))+
  geom_point()+
  geom_line()+
  facet_grid(itn_cov ~ ivm_cov)

#ggsave(Sx_plot, file = "plots/Sxtot_ivm_llin.svg")
#ggsave(Ex_plot, file = "plots/Extot_ivm_llin.svg")
#ggsave(Ix_plot, file = "plots/Ixtot_ivm_llin.svg")
#
out.df2$mv
mvx_prop_high_endo <- ggplot(out.df2, aes(x = ivm_cov, y = mvxtot/mv, fill = as.factor(itn_cov)))+
  geom_bar(position = position_dodge(), stat="identity")+
  ggtitle("High endophagy")+
  ylim(0, 1)

mvx_prop_low_endo <- ggplot(out.df3, aes(x = ivm_cov, y = mvxtot/mv, fill = as.factor(itn_cov)))+
  geom_bar(position = position_dodge(), stat="identity")+
  ggtitle("Low endophagy")+
  ylim(0, 1)

plot_grid(mvx_prop_high_endo, mvx_prop_low_endo)

#ggsave(mvx_plot, file = "plots/mvxtot_ivm_llin.svg")

out.df2.summary <- out.df2 %>%
  filter(itn_cov != 1 & ivm_cov != 1) %>%
  group_by(itn_cov, ivm_cov) %>%
  summarise(Sx_dead_tot = sum(Sx_dead),
            Ex_dead_tot = sum(Ex_dead),
            Ix_dead_tot = sum(Ix_dead),
            mvx_dead_tot = sum(mvx_dead))

unique(out.df2.summary$ivm_cov)

mvx_dead_plot <- ggplot(out.df2.summary, aes(x = ivm_cov, y = mvx_dead_tot, fill = as.factor(itn_cov)))+
  geom_bar(position = position_dodge(), stat="identity")+
  #geom_point()+
  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed mosquitoes",
       fill = "LLIN coverage")+
  theme_minimal()+
  ylim(0, 2e5)


#Sx_dead_plot <- ggplot(out.df2.summary, aes(x = ivm_cov, y = Sx_dead_tot, fill = as.factor(itn_cov)))+
#  geom_bar(stat = "identity")+
#  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed Sx mosquitoes",
#       fill = "LLIN coverage")+
#  theme_minimal()+
#  ylim(0, 6e5)+
#  ggtitle("High endophagy")
#
#Ex_dead_plot <- ggplot(out.df2.summary, aes(x = ivm_cov, y = Ex_dead_tot, fill = as.factor(itn_cov)))+
#  geom_bar(stat = "identity")+
#  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed Ex mosquitoes",
#       fill = "LLIN coverage")+
#  theme_minimal()+
#  ylim(0, 4e4)

Ix_dead_plot <- ggplot(out.df2.summary, aes(x = ivm_cov, y = Ex_dead_tot,  fill = as.factor(itn_cov)))+
  geom_bar(position="dodge", stat="identity")+
  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed Ix mosquitoes",
       fill = "LLIN coverage")+
  theme_minimal()+
  ylim(0, 1e4)+
  ggtitle("High endophagy")

plot_grid(Ix_dead_plot, mvx_dead_plot, ncol = 2)



#plots with low endo
out.df3.summary <- out.df3 %>%
  filter(itn_cov != 1 & ivm_cov != 1) %>%
  group_by(itn_cov, ivm_cov) %>%
  summarise(Sx_dead_tot = sum(Sx_dead),
            Ex_dead_tot = sum(Ex_dead),
            Ix_dead_tot = sum(Ix_dead),
            mvx_dead_tot = sum(mvx_dead))

mvx_dead_plot_low_endo <- ggplot(out.df3.summary, aes(x = ivm_cov, y = mvx_dead_tot, fill = as.factor(itn_cov)))+
  geom_bar(position="dodge", stat="identity")+
  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed mosquitoes",
       fill = "LLIN coverage")+
  theme_minimal()+
  ylim(0, 2e5)

#Sx_dead_plot_low_endo <- ggplot(out.df3.summary, aes(x = ivm_cov, y = Sx_dead_tot, fill = as.factor(itn_cov)))+
#  geom_bar(stat = "identity")+
#  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed Sx mosquitoes",
#       fill = "LLIN coverage")+
#  theme_minimal()+
#  ylim(0, 6e5)+
#  ggtitle("Low endophagy")
#
#Ex_dead_plot_low_endo <- ggplot(out.df3.summary, aes(x = ivm_cov, y = Ex_dead_tot, fill = as.factor(itn_cov)))+
#  geom_bar(stat = "identity")+
#  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed Ex mosquitoes",
#       fill = "LLIN coverage")+
#  theme_minimal()+
#  ylim(0, 4e4)

Ix_dead_plot_low_endo <- ggplot(out.df3.summary, aes(x = ivm_cov, y = Ex_dead_tot, fill = as.factor(itn_cov)))+
  geom_bar(position="dodge", stat="identity")+
  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed Ix mosquitoes",
       fill = "LLIN coverage")+
  theme_minimal()+
  ylim(0, 1e04)+
  ggtitle("Low endophagy")

#plot_grid(Sx_dead_plot_low_endo, Ex_dead_plot_low_endo, Ix_dead_plot_low_endo, mvx_dead_plot_low_endo, ncol = 2,
#          nrow = 2)

dead_mosq_endo <- plot_grid(
          Ix_dead_plot, Ix_dead_plot_low_endo,
          mvx_dead_plot, mvx_dead_plot_low_endo,
          ncol = 2, nrow = 2)

#ggsave(dead_mosq_endo, file = "plots/dead_mosq_endo.svg")


#look at EIR: absolute change and relative

#absolute EIR

#EIR over time in each model

#just checking dynamics with no ivermectin
out.df2 %>%
  filter(ivm_cov ==  0) %>%
  ggplot(aes(x = t, y = EIR_tot, col = as.factor(itn_cov)))+
  geom_line()+
  ylim(0, 300)


#have a higher EIR when you introduce bed nets because more transmission is left over
#when bites_indoors and bites_bed are low
out.df3 %>%
  filter(ivm_cov ==  0) %>%
  ggplot(aes(x = t, y = EIR_tot, col = as.factor(itn_cov)))+
  geom_line()+
  ylim(0, 300)

#compute mean EIR over the time period (10 years)

#introducing ivermectin on d3100.

no_ivm <- out.df2$ivm_cov==0

out.df2.eir <- out.df2 %>%
  mutate(ivm_on = case_when(2735 <= t & t <= 2795 + 23 ~ "0", #corr time pre ivm distrib
                            3100 <= t & t <= 3160 + 23 ~ "1", #ivm distrib
                             TRUE ~ "3")) %>% #after ivm distrib
  filter(itn_cov != 0 & ivm_cov != 0) %>%
  group_by(itn_cov, ivm_cov, ivm_on) %>%
  summarise(mean_EIR = mean(EIR_tot)) %>%
  mutate(endo = "high_endo")

out.df3.eir <- out.df3 %>%
  mutate(ivm_on = case_when(t < 3100 ~ "0", #before ivm distrib
                             TRUE ~ "1")) %>% #after ivm distrib
  filter(itn_cov != 0 & ivm_cov != 0) %>%
  group_by(itn_cov, ivm_cov, ivm_on) %>%
  summarise(mean_EIR = mean(EIR_tot)) %>%
  mutate(endo = "low_endo")


out.df.eir.endo <- rbind(out.df2.eir, out.df3.eir)

ggplot(out.df2.eir, aes(x  = ivm_on, y = mean_EIR, col = as.factor(ivm_cov)))+
  geom_point()+
  #geom_line(aes(col = as.factor(ivm_cov)), group = 1)+
  ylim(0, 150)+
  ylab("Mean EIR")+
  xlab("Ivermectin off/on")+
  ggtitle("High endophagy")+
  facet_wrap(vars(itn_cov))



high_endo_eir <- ggplot(out.df3.eir, aes(x  = ivm_cov, y = mean_EIR, fill = as.factor(itn_cov)))+
  geom_bar(position = "dodge", stat = "identity")+
  ylim(0, 150)+
  ylab("Mean EIR across 10y model run")+
  xlab("Ivermectin coverage")+
  labs(fill = "LLIN coverage")+
  ggtitle("Low endophagy")+
  theme_minimal()

low_endo_eir <- ggplot(out.df3.eir, aes(x  = ivm_cov, y = mean_EIR, fill = as.factor(itn_cov)))+
  geom_bar(position = "dodge", stat = "identity")+
  ylim(0, 150)+
  ylab("Mean EIR across 10y model run")+
  xlab("Ivermectin coverage")+
  labs(fill = "LLIN coverage")+
  ggtitle("Low endophagy")+
  theme_minimal()

cowplot::plot_grid(high_endo_eir, low_endo_eir)

ggplot(out.df2, aes(x = ivm_cov, y = EIR_tot, fill = as.factor(itn_cov)))+
  geom_bar(position = position_dodge(), stat="identity")+
  #geom_point()+
  labs(x = "Ivermectin coverage", y = "Number of \n ivm-killed mosquitoes",
       fill = "LLIN coverage")+
  theme_minimal()

feeding_ivm_uptake <- cowplot::plot_grid(mvx_plot, mvx_plot_low_endo, ncol = 2)
ggsave(feeding_ivm_uptake, file = "plots/endophagy_ivm_llin.svg")
ggsave(feeding_ivm_uptake, file = "plots/endophagy_ivm_llin.pdf",
       device = cairo_pdf,
       width = 297,
       height = 210,
       units = "mm")



prop_fed_ivm_plot <- out.df2 %>%
  mutate(prop_fed_ivm = mvxtot/mv) %>%
  ggplot()+
  aes(x = t, y = prop_fed_ivm) +
  geom_point()+
  geom_line()+
  facet_grid(itn_cov ~ ivm_cov)
ggsave(prop_fed_ivm_plot, file = "plots/prop_fed_ivm_plot.svg")

out.df2 %>%
  ggplot()+
  aes(x = t, y = mv)+
  geom_point()+
  geom_line()+
  facet_grid(itn_cov ~ ivm_cov)


range(out.df2$Sxtot)
range(out.df2$Extot)
range(out.df2$Ixtot)

ggplot(out.df2, aes(x = t, y = Sxtot))+
  geom_point()

require(cowplot)

plot1 <-ggplot(out.df2, aes(x = t, y = avhc, col= as.factor(itn_cov)))+
  geom_point()+
  ylim(0, 1)

plot2 <- ggplot(out.df2, aes(x = t, y = FOIv, col = as.factor(itn_cov)))+
  geom_point()



avhc_foi_plot <- plot_grid(plot1, plot2)
ggsave(avhc_foi_plot, file = "plots/avhc_foi_plot.svg")


#looking at relationship between avhc and r_ITN0 and d_ITN0. Default bites_indoors and bites_bed parameters
out.df4 <- read.csv("data/out_df_4.csv", header = TRUE)

a <- ggplot(out.df4, aes(x = t, y = avhc, col = as.factor(d_ITN0)))+
  geom_line()

b <- ggplot(out.df4, aes(x = t, y = avhc, col = as.factor(r_ITN0)))+
  geom_line()

cowplot::plot_grid(a, b)

ggplot(out.df4, aes(x = t, y = mv, col = as.factor(d_ITN0)))+
  geom_line()

ggplot(out.df4, aes(x = t, y = d_ITN, col = as.factor(d_ITN0)))+
  geom_line()+
  ylim(0, 1)

ggplot(out.df4, aes(x = t, y = r_ITN, col = as.factor(d_ITN0)))+
  geom_line()

ggplot(out.df4, aes(x = t, y = s_ITN, col = as.factor(d_ITN0)))+
  geom_line()

out.df4.summary <- out.df4 %>%
  group_by(d_ITN0, r_ITN0, ivm_cov) %>%
  summarise(Sx_dead_tot = sum(Sx_dead),
            Ex_dead_tot = sum(Ex_dead),
            Ix_dead_tot = sum(Ix_dead),
            mvx_dead_tot = sum(mvx_dead))

ggplot(out.df4.summary, aes(x = ivm_cov, y = mvx_dead_tot), col = as.factor(d_ITN0))+
  geom_bar(stat = "identity")

#this doesn't make sense to me - surely s_ITN should vary depending on the starting d_ITN0
ggplot(out.df4, aes(x = t, y = s_ITN, col = as.factor(d_ITN0)))+
  geom_line()



ggplot(out.df4, aes(x = t, y = Q, col = as.factor(d_ITN0)))+
  geom_line()

out.df4 %>%
  select(t, d_ITN0, r_ITN0, s_ITN, d_ITN, r_ITN, ITN_decay, itn_loss)

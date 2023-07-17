#plots for the llin_ivm.R script
require(tidyverse)

#av_mosq and avhc under different LLIN coverages
out.df <- read.csv("data/out_df.csv", header = TRUE)

#PLOT 1
#plot of av_mosq_sum (across the intervention categories changing with ITN cov). Difficult to interpret because summing
ggplot(out.df, aes(x = t, y = av_mosq_sum, col = as.factor(itn_cov)))+
  geom_line(linewidth = 1)

#PLOT 2
#plot of avhc over time - this is the normalised value and what clearly affects the ivm uptake
ggplot(out.df, aes(x = t, y = avhc, col = as.factor(itn_cov)))+
  geom_line(linewidth = 1)


#looking at how uptake of ivermectin by mosquitoes changes with LLIN cov and IVM cov
out.df2 <- read.csv("data/out_df_2.csv", header = TRUE)
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

ggsave(Sx_plot, file = "plots/Sxtot_ivm_llin.svg")
ggsave(Ex_plot, file = "plots/Extot_ivm_llin.svg")
ggsave(Ix_plot, file = "plots/Ixtot_ivm_llin.svg")

mvx_plot <- ggplot(out.df2, aes(x = t, y = mvxtot, col = as.factor(is_neg)))+
  geom_point()+
  geom_line()+
  facet_grid(itn_cov ~ ivm_cov)
ggsave(mvx_plot, file = "plots/mvxtot_ivm_llin.svg")

prop_fed_ivm_plot <- out.df2 %>%
  mutate(prop_fed_ivm = mvxtot/mv) %>%
  ggplot()+
  aes(x = t, y = prop_fed_ivm) +
  geom_point()+
  geom_line()+
  facet_grid(itn_cov ~ ivm_cov)
ggsave(prop_fed_ivm_plot, file = "plots/prop_fed_ivm_plot.svg")


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


plot3 <- ggplot(out.df2, aes(x = t, y = lag_FOIv, col = as.factor(itn_cov)))+
  geom_point()

avhc_foi_plot <- plot_grid(plot1, plot2, plot3)

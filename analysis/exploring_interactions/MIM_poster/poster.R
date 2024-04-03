

require(tidyverse)

antag_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/antag.rds")
add_mod <- readRDS(file = "analysis/exploring_interactions/MIM_poster/add.rds")

mods <- rbind(antag_mod, add_mod)

res <- unique(antag_mod$d_ITN0)
itn_cov <- unique(antag_mod$itn_cov)

mods_dynamics <- mods %>%
  filter(d_ITN0 == res[1] & model %in% c("antag_LLIN", "antag_LLIN_IVM") & itn_cov == 0.8)

net_seq <- seq(100, 3650, by = 3*365)
IVM_begin1 <- net_seq[3]+(6*30) # 6 months into new net distribution
mda_int <- 30
IVM_start1 <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)

lines <- c("antag_LLIN" = "dotted", "antag_LLIN_IVM" = "solid")

sp_pals <- c('#1b9e77','#d95f02','#7570b3','#e7298a')

ggplot(mods_dynamics, aes(x = t/365, y = slide_prev0to5*100, linetype = as.factor(model), col = as.factor(species)))+
  geom_line(linewidth = 1)+
  theme_minimal()+
  labs(y = "Slide prevalence in under 5s (%)",
       x = "Time (years)")+
  scale_linetype_manual(values = lines, name = "Intervention", labels = c("LLIN", "LLIN & endectocide"))+
  ylim(0, 48)+
  annotate("segment", x = IVM_start1[1]/365, xend = IVM_start1[1]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[2]/365, xend = IVM_start1[2]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = IVM_start1[3]/365, xend = IVM_start1[3]/365, y = 35, yend = 31, colour = "black", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[1]/365, xend = net_seq[1]/365, y = 45, yend = 41, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[2]/365, xend = net_seq[2]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[3]/365, xend = net_seq[3]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate("segment", x = net_seq[4]/365, xend = net_seq[4]/365, y = 35, yend = 31, colour = "#1f78b4", arrow = arrow(length = unit(0.01, "npc")))+
  annotate(geom = "text", x = (net_seq[1]/365)+1, y = 47, label = "LLIN distribution (every 3 years)", col = "#1f78b4")+
  annotate(geom = "text", x = (IVM_start1[1]/365)+(420/365), y = 36, label = "Endectocide MDA campaign")+
  scale_colour_manual(name = "Anopheles species", labels = c("arabiensis",
                                                             "funestus",
                                                             "gambiae",
                                                             "stephensi"),
                      values = sp_pals)+
  theme(legend.position = c(0.35, 0.15), legend.direction = "horizontal")


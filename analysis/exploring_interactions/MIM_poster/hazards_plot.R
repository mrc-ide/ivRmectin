require(tidyverse)
ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")

ivm_haz <- as.data.frame(ivm_haz)

hr_plot <- ggplot(ivm_haz, aes(x = Day, y = IVM_300_3_HS))+
  geom_point(size = 2)+
  labs(x = expression("Day of bloodmeal after first dose of ivermectin-like drug (3x300 " * mu * "g/kg)"),
       y = "Hazard ratio")+
  theme_minimal()+
  geom_vline(xintercept = 23, lty = "dashed", col = "red")+
  ylim(0, 10)+
  theme(text = element_text(size = 20))

ggsave(hr_plot, file = "analysis/exploring_interactions/MIM_poster/hazards.svg",
       width = 30.31,
       height = 14.22,
       units = "cm")

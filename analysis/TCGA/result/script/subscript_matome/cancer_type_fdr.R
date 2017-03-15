library(dplyr)
library(ggplot2)

source("../../../conf/plot_config.R")
  
cacner_type_fdr <- read.table("../temporary/cancer_type_fdr.txt", sep = "\t", header = TRUE)

ggplot(cacner_type_fdr, aes(x = Cancer_Type, y = FDR)) +
  geom_bar(stat = "identity", fill = "#6baed6") +
  my_theme() +
  labs(x = "Cancer type") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0))


ggsave("../figure/cancer_type_fdr.tiff", width = 10, height = 5, dpi = 600, units = "cm")


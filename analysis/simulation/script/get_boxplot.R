#! /usr/local/package/r/3.2.5/bin/R 

library(dplyr)
library(ggplot2)

source("../../conf/plot_config.R")

sim1 <- read.table("../output/sim1.output.txt", sep = "\t", header = TRUE)
sim1$Splicing_Num <- factor(sim1$Splicing_Num)

ggplot(sim1, aes(x = Splicing_Num, y = BF, fill = Is_Active)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.3) +
  my_theme() +
  labs(x = "#splicing pattern", y = "log(Bayes Factor)", fill = "") +
  theme(legend.position = "bottom") 

ggsave("../output/sim1.boxplot.tiff", width = 7, height = 7, dpi = 600, units = "cm")

sim2 <- read.table("../output/sim2.output.txt", sep = "\t", header = TRUE)
sim2$Splicing_Num <- factor(sim2$Splicing_Num)

ggplot(sim2, aes(x = Splicing_Num, y = BF, fill = Is_Active)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.3) +
  my_theme() +
  labs(x = "#mutation", y = "log(Bayes Factor)", fill = "") +
  theme(legend.position = "bottom")

ggsave("../output/sim2.boxplot.tiff", width = 7, height = 7, dpi = 600, units = "cm")


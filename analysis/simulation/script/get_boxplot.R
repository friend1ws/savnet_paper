#! /usr/local/package/r/3.2.5/bin/R 

library(dplyr)
library(ggplot2)

sim1 <- read.table("../output/sim1.output.txt", sep = "\t", header = TRUE)
sim1$Splicing_Num <- factor(sim1$Splicing_Num)

ggplot(sim1, aes(x = Splicing_Num, y = BF, fill = Is_Active)) + 
  geom_boxplot() +
  theme_minimal() +
  labs(x = "#splicing", y = "log(Bayes Factor)") +
  theme(legend.position = "bottom")

ggsave("../output/sim1.boxplot.pdf", width = 5, height = 5)

sim2 <- read.table("../output/sim2.output.txt", sep = "\t", header = TRUE)
sim2$Splicing_Num <- factor(sim2$Splicing_Num)

ggplot(sim2, aes(x = Splicing_Num, y = BF, fill = Is_Active)) + 
  geom_boxplot() +
  theme_minimal() +
  labs(x = "#mutation", y = "log(Bayes Factor)") +
  theme(legend.position = "bottom")

ggsave("../output/sim2.boxplot.pdf", width = 5, height = 5)


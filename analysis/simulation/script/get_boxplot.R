#! /usr/local/package/r/3.2.5/bin/R 

library(dplyr)
library(ggplot2)

source("../../conf/plot_config.R")

sim1 <- read.table("../output/sim1.output.txt", sep = "\t", header = TRUE)
sim1$Is_Active <- factor(sim1$Is_Active, levels = c("active", "inactive"), labels = c("True association", "False association"))

sim1$Splicing_Num <- factor(sim1$Splicing_Num)

sim1_sensitivity <- sim1 %>% 
    mutate(Is_Call = ifelse(BF >= 3.0, 1, 0)) %>% 
    group_by(Is_Active, Splicing_Num) %>% 
    summarize(Sensitivity = mean(Is_Call))


ggplot() + 
  geom_boxplot(data = sim1, aes(x = Splicing_Num, y = BF, fill = Is_Active), size = 0.3, outlier.size = 0.3) +
  geom_point(data = sim1_sensitivity %>% filter(Is_Active == "True association"), 
             aes(x = as.numeric(Splicing_Num), y = 50 * (Sensitivity - 0.2)), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  geom_line(data = sim1_sensitivity %>% filter(Is_Active == "True association"), 
            aes(x = as.numeric(Splicing_Num), y = 50 * (Sensitivity - 0.2)), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  geom_abline(intercept = 3.0, slope = 0, colour = "#d73027", alpha = 0.6, linetype = "longdash") + 
  my_theme() +
  labs(x = "Number of possible association", y = "Log10(Bayes Factor)", fill = "") +
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = c("True association" = "#b2182b", "False association" = "#bababa")) +
  scale_y_continuous(limits = c(-10, 40), breaks = c(-10, 0, 3, 10, 20, 30, 40), sec.axis = sec_axis(~ . * (1 / 50) + 0.2, name = "Sensitivity")) +


ggsave("../output/sim1.boxplot.tiff", width = 8, height = 6, dpi = 600, units = "cm")



sim2 <- read.table("../output/sim2.output.txt", sep = "\t", header = TRUE)
sim2$Is_Active <- factor(sim2$Is_Active, levels = c("active", "inactive"), labels = c("True association", "False association"))

sim2$Splicing_Num <- factor(sim2$Splicing_Num)

sim2_sensitivity <- sim2 %>%
    mutate(Is_Call = ifelse(BF >= 3.0, 1, 0)) %>%
    group_by(Is_Active, Splicing_Num) %>%
    summarize(Sensitivity = mean(Is_Call))


ggplot() + 
  geom_boxplot(data = sim2, aes(x = Splicing_Num, y = BF, fill = Is_Active), size = 0.3, outlier.size = 0.3) +
  geom_point(data = sim2_sensitivity %>% filter(Is_Active == "True association"),
             aes(x = as.numeric(Splicing_Num), y = 50 * (Sensitivity - 0.2)), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  geom_line(data = sim2_sensitivity %>% filter(Is_Active == "True association"),
            aes(x = as.numeric(Splicing_Num), y = 50 * (Sensitivity - 0.2)), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  geom_abline(intercept = 3.0, slope = 0, colour = "#d73027", alpha = 0.6, linetype = "longdash") +
  my_theme() +
  labs(x = "Number of associated variants", y = "Log10(Bayes Factor)", fill = "") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("True association" = "#b2182b", "False association" = "#bababa")) +
  scale_y_continuous(limits = c(-10, 40), breaks = c(-10, 0, 3, 10, 20, 30, 40), sec.axis = sec_axis(~ . * (1 / 50) + 0.2, name = "Sensitivity")) +

ggsave("../output/sim2.boxplot.tiff", width = 8, height = 6, dpi = 600, units = "cm")



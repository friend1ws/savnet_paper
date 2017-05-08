library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

source("../../../conf/plot_config.R")

pos_fdr <- read.table("../temporary/position_fdr.txt", sep = "\t", header = TRUE)


total_pos_fdr <- pos_fdr %>% 
  group_by(Is_Original, Mutation_Type, Rel_Pos) %>% 
  summarize(total_count = sum(Count)) %>%
  spread(key = Is_Original, value = total_count) %>%
  mutate(pfdr = permutation / (100 * original))

total_pos_fdr[total_pos_fdr$pfdr > 1, "pfdr"] <- 1.0
# total_pos_fdr[total_pos_fdr$Is_Original == "permutation","total_count"] <- 
#   total_pos_fdr$total_count[total_pos_fdr$Is_Original == "permutation"] / 100

##########
total_pos_fdr_donor <- total_pos_fdr %>% 
  filter(Mutation_Type == "splicing donor disruption")

total_pos_fdr_donor$Rel_Pos2 <- 
  factor(total_pos_fdr_donor$Rel_Pos, 
         labels = c("-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5",
                    "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15"))


##########

##########
total_pos_fdr_acceptor <- total_pos_fdr %>% 
  filter(Mutation_Type == "splicing acceptor disruption")

total_pos_fdr_acceptor$Rel_Pos2 <- 
  factor(total_pos_fdr_acceptor$Rel_Pos, 
         labels = rev(c("-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5",
                        "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15")))


##########

df_donor <- pos_fdr %>% 
  filter(Mutation_Type == "splicing donor disruption") %>% 
  group_by(Is_Original, Rel_Pos) %>% 
  summarize(Total_Count = sum(Count)) 

df_donor[df_donor$Is_Original == "permutation", "Total_Count"] <- 
  df_donor[df_donor$Is_Original == "permutation", "Total_Count"] / 100

df_donor$Rel_Pos2 <- 
  factor(df_donor$Rel_Pos, 
         labels = c("-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5",
                    "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15"))

df_donor$Is_Original2 <- 
  factor(df_donor$Is_Original, labels = c("Called SAV count", "Estimated false positive count"))

df_region <- data.frame(
  xmin = c(0, 5.5),
  xmax = c(5.5, 21),
  ymin = c(-Inf, -Inf),
  ymax = c(Inf, Inf),
  region = c("exon", "intron")
)

g_donor_2 <- ggplot() +
  geom_bar(data = df_donor, aes(x = Rel_Pos2, y = Total_Count, fill = Is_Original2),
           stat = "identity", position = "dodge") +
  geom_point(data = total_pos_fdr_donor, aes(x = Rel_Pos2, y = 5000 * pfdr), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  geom_line(data = total_pos_fdr_donor, aes(x = Rel_Pos, y = 5000 * pfdr), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  ggtitle("Donor disruption") +
  labs(x = "Position", y = "Number") +
  geom_vline(xintercept = 5.5, colour="#d73027", linetype = "longdash") +
  my_theme() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"),
        panel.grid.major.x = element_blank()) +
  scale_fill_manual(values = c("Called SAV count" = "#66bd63", "Estimated false positive count" = "#bf812d")) +
  scale_x_discrete(labels = 
                     c("-5", rep("", 3), "-1 ", " +1", rep("", 3), "+5", rep("", 4), "+10", rep("", 4), "+15")) +
  scale_y_continuous(limits = c(0, 4000), sec.axis = sec_axis(~ . * (1 / 5000), name = "")) +
  guides(fill = FALSE)

##########

df_acceptor <- pos_fdr %>% 
  filter(Mutation_Type == "splicing acceptor disruption") %>% 
  group_by(Is_Original, Rel_Pos) %>% 
  summarize(Total_Count = sum(Count)) 

df_acceptor[df_acceptor$Is_Original == "permutation", "Total_Count"] <- 
  df_acceptor[df_acceptor$Is_Original == "permutation", "Total_Count"] / 100

df_acceptor$Rel_Pos2 <- 
  factor(df_acceptor$Rel_Pos, 
         labels = rev(c("-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5",
                        "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15")))

df_acceptor$Is_Original2 <- 
  factor(df_acceptor$Is_Original, labels = c("Called SAV count", "Estimated false positive count"))

g_acceptor_2 <- ggplot() + 
  geom_bar(data = df_acceptor, aes(x = Rel_Pos2, y = Total_Count, fill = Is_Original2),
           stat = "identity", position = "dodge") +
  geom_point(data = total_pos_fdr_acceptor, aes(x = Rel_Pos2, y = 5000 * pfdr), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  geom_line(data = total_pos_fdr_acceptor, aes(x = Rel_Pos, y = 5000 * pfdr), size = 0.6, alpha = 0.8, colour = "#7570b3") +
  geom_vline(xintercept = 15.5, colour="#d73027", linetype = "longdash") +
  ggtitle("Acceptor disruption") +
  labs(x = "Position", y = "") +
  my_theme() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"),
        panel.grid.major.x = element_blank()) +
  scale_fill_manual(values = c("Called SAV count" = "#66bd63", "Estimated false positive count" = "#bf812d")) +
  scale_x_discrete(labels = 
                     c("+15", rep("", 4), "+10", rep("", 4), "+5", 
                       rep("", 3), "+1 ", " -1", rep("", 3), "-5")) +
  scale_y_continuous(limits = c(0, 4000), sec.axis = sec_axis(~ . * (1 / 5000), name = "Position-wise FDR")) +
  guides(fill = FALSE)



g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

g_dummy_for_legend <- 
  ggplot(df_acceptor,
         aes(x = Rel_Pos2, y = Total_Count, fill = Is_Original2)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Called SAV count" = "#66bd63", "Estimated false positive count" = "#bf812d")) +
  my_theme() + 
  theme(legend.position = "bottom") +
  labs(fill = "")


g_donor_acceptor_2 <- plot_grid(g_donor_2, g_acceptor_2, ncol = 2, align = "h")

plot_grid(g_donor_acceptor_2, ggdraw() + draw_label(" ", size = 7), g_legend(g_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.10, 0.10))

ggsave("../figure/position_fpnum_fdr.tiff", width = 12, height = 5.0, dpi = 600, unit = "cm")




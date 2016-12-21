library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

pos_fdr <- read.table("../matome/position_fdr.txt", sep = "\t", header = TRUE)

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

g_donor <- ggplot(total_pos_fdr_donor,
                  aes(x = Rel_Pos2, y = pfdr)) +
  ggtitle("Donor") +
  labs(x = "Intron Position", y = "Position-Wise FDR") +
  geom_rect(xmin = 0, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#377eb8") +
  geom_rect(xmin = 5.5, xmax = 21, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#e41a1c") +
  geom_bar(fill = "#984ea3", stat = "identity") +
  ylim(c(0, 0.8)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(labels = 
                     c("-5", rep("", 3), "-1", "+1", "+2", rep("", 2), "+5", rep("", 4), "+10", rep("", 4), "+15"))


##########

##########
total_pos_fdr_acceptor <- total_pos_fdr %>% 
  filter(Mutation_Type == "splicing acceptor disruption")

total_pos_fdr_acceptor$Rel_Pos2 <- 
  factor(total_pos_fdr_acceptor$Rel_Pos, 
         labels = rev(c("-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5",
                        "+6", "+7", "+8", "+9", "+10", "+11", "+12", "+13", "+14", "+15")))

g_acceptor <- ggplot(total_pos_fdr_acceptor,
                     aes(x = Rel_Pos2, y = pfdr)) +
  ggtitle("Acceptor") +
  labs(x = "Intron Position", y = "Position-Wise FDR") +
  geom_rect(xmin = 0, xmax = 15.5, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#e41a1c") +
  geom_rect(xmin = 15.5, xmax = 21, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#377eb8") +
  geom_bar(fill = "#984ea3", stat = "identity") +
  ylim(c(0, 0.8)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(labels = 
                     c("+15", rep("", 4), "+10", rep("", 4), "+5", 
                       rep("", 2), "+2", "+1", "-1", rep("", 3), "-5"))


##########

plot_grid(g_donor, g_acceptor, ncol = 2, align = "h")


ggsave("../matome/position_fdr.png", width = 8, height = 3, unit = "in")

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
  factor(df_donor$Is_Original, labels = c("#Called", "#Estimated False Positives"))

g_donor_2 <- ggplot(df_donor,
                    aes(x = Rel_Pos2, y = Total_Count, fill = Is_Original2)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Donor") +
  labs(x = "Intron Position", y = "Count") +
  geom_rect(xmin = 0, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#377eb8") +
  geom_rect(xmin = 5.5, xmax = 21, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#e41a1c") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = 
                     c("-5", rep("", 3), "-1", "+1", "+2", rep("", 2), "+5", rep("", 4), "+10", rep("", 4), "+15")) +
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
  factor(df_acceptor$Is_Original, labels = c("#Called", "#Estimated False Positives"))

g_acceptor_2 <-
  ggplot(df_acceptor,
         aes(x = Rel_Pos2, y = Total_Count, fill = Is_Original2)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Acceptor") +
  labs(x = "Intron Position", y = "Count") +
  geom_rect(xmin = 0, xmax = 15.5, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#e41a1c") +
  geom_rect(xmin = 15.5, xmax = 21, ymin = -Inf, ymax = Inf, alpha = 0.002, fill = "#377eb8") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = 
                     c("+15", rep("", 4), "+10", rep("", 4), "+5", 
                       rep("", 2), "+2", "+1", "-1", rep("", 3), "-5")) +
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
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "bottom") +
  labs(fill = "")


g_donor_acceptor_2 <- plot_grid(g_donor_2, g_acceptor_2, ncol = 2, align = "h")

plot_grid(g_donor_acceptor_2, g_legend(g_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.2))

ggsave("../matome/position_fpnum.png", width = 8, height = 3, unit = "in")



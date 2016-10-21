library(dplyr)
library(ggplot2)


gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
base_col <- c("A" = gg_color_hue6[3], "C" = gg_color_hue6[5], "G" = gg_color_hue6[2], "T" = gg_color_hue6[1])



target_motif_info <- read.table("../matome/omega.motif_summary.txt", header = TRUE, sep = "\t")

target_motif_count <- target_motif_info %>% 
  group_by(Rel_Pos, Alt_Base, Motif_Type) %>% 
  summarize(count = n())

##########
# splicing donor motif (target)
pos_colour <- rep("grey30", 8)
pos_colour[3:4] <- "red"

ggplot(target_motif_count %>% filter(Motif_Type == "splicing donor disruption"),
       aes(x = Rel_Pos, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "#SNV", fill = "alternative base") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = pos_colour, size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:8, 
                   labels = c("A", "G", "G", "T", "R", "A", "G", "T")) +
  scale_fill_manual(values = base_col)
  
ggsave("../matome/donor_dis_target.png", width = 10, height = 6)

##########
# splicing acceptor motif

pos_colour <- rep("grey30", 8)
pos_colour[7:8] <- "red"

ggplot(target_motif_count %>% filter(Motif_Type == "splicing acceptor disruption"),
       aes(x = Rel_Pos, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "#SNV", fill = "alternative base") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = pos_colour, size = rel(2)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:9, 
                   labels = c("Y", "Y", "Y", "Y", "N", "C", "A", "G", "G")) +
  scale_fill_manual(values = base_col)

ggsave("../matome/acceptor_dis_target.png", width = 10, height = 6)




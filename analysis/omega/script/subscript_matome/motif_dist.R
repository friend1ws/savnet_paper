library(dplyr)
library(ggplot2)
library(cowplot)

gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
base_col <- c("A" = gg_color_hue6[3], "C" = gg_color_hue6[5], "G" = gg_color_hue6[2], "T" = gg_color_hue6[1])


target_motif_info <- read.table("omega.motif_summary.txt", header = TRUE, sep = "\t")
# target_motif_info <- read.table("../matome/omega.motif_summary.txt", header = TRUE, sep = "\t")

target_motif_count <- target_motif_info %>% 
  group_by(Rel_Pos, Ref_Base, Alt_Base, Motif_Type) %>% 
  summarize(count = n())


##########





##########
# splicing donor motif (target)
# pos_colour <- rep("grey30", 8)
# pos_colour[3:4] <- "red"

target_motif_count_dd <- target_motif_count %>% filter(Motif_Type == "splicing donor disruption")

target_motif_count_dd$Rel_Pos2 <- rep(NA, nrow(target_motif_count_dd))
target_motif_count_dd$Rel_Pos2[target_motif_count_dd$Rel_Pos <= 2] <- 
  target_motif_count_dd$Rel_Pos[target_motif_count_dd$Rel_Pos <= 2] - 3

target_motif_count_dd$Rel_Pos2[target_motif_count_dd$Rel_Pos >= 3] <- 
  target_motif_count_dd$Rel_Pos[target_motif_count_dd$Rel_Pos >= 3] - 2


target_motif_count_dd$Rel_Pos2 <- 
  factor(target_motif_count_dd$Rel_Pos2, levels = c("-2", "-1", "1", "2", "3", "4", "5", "6"))

pos_colour <- rep("grey30", 8)
pos_colour[3:4] <- "red"

p_dd <- ggplot(target_motif_count_dd,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "reference base", y = "#SNV", fill = "alternative base") +
  ylim(c(0, 2250)) +
  theme_minimal() +
  ggtitle("splicing donor disruption") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  guides(fill = FALSE)

# ggsave("../matome/donor_dis_target.png", width = 10, height = 6)

##########
# splicing acceptor motif


target_motif_count_ad <- target_motif_count %>% filter(Motif_Type == "splicing acceptor disruption")

target_motif_count_ad$Rel_Pos2 <- rep(NA, nrow(target_motif_count_ad))
target_motif_count_ad$Rel_Pos2[target_motif_count_ad$Rel_Pos <= 8] <- 
  9 - target_motif_count_ad$Rel_Pos[target_motif_count_ad$Rel_Pos <= 8] 

target_motif_count_ad$Rel_Pos2[target_motif_count_ad$Rel_Pos >= 9] <- -1


target_motif_count_ad$Rel_Pos2 <- 
  factor(target_motif_count_ad$Rel_Pos2, levels = c("8", "7", "6", "5", "4", "3", "2", "1", "-1"))


p_ad <- ggplot(target_motif_count_ad,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Rel_Pos2, nrow = 1) +
  labs(x = "reference base", y = "#SNV", fill = "alternative base") +
  ylim(c(0, 2250)) +
  theme_minimal() +
  ggtitle("splicing acceptor disruption") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
scale_fill_manual(values = base_col) +
  guides(fill = FALSE)


##########
# legend

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p_dummy_for_legend <- ggplot(target_motif_count_dd,
                             aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "reference base", y = "#SNV", fill = "alternative base") +
  theme_minimal() +
  ggtitle("splicing donor disruption") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col)

# ggsave("../matome/acceptor_dis_target.png", width = 10, height = 6)

# p_dd_ad <- plot_grid(p_dd, p_ad, ncol = 2, align = "h", rel_widths = c(0.9, 1))

# plot_grid(p_dd_ad, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.1))


# ggsave("donor_acceptor_disruption.png", width = 10, height = 3.5)


##########
# splicing donor creation
target_motif_count_dc <- target_motif_count %>% filter(Motif_Type == "splicing donor creation")

target_motif_count_dc$Rel_Pos2 <- rep(NA, nrow(target_motif_count_dc))
target_motif_count_dc$Rel_Pos2[target_motif_count_dc$Rel_Pos <= 2] <- 
  target_motif_count_dc$Rel_Pos[target_motif_count_dc$Rel_Pos <= 2] - 3

target_motif_count_dc$Rel_Pos2[target_motif_count_dc$Rel_Pos >= 3] <- 
  target_motif_count_dc$Rel_Pos[target_motif_count_dc$Rel_Pos >= 3] - 2

target_motif_count_dc$Rel_Pos2 <- 
  factor(target_motif_count_dc$Rel_Pos2, levels = c("-2", "-1", "1", "2", "3", "4", "5", "6"))

target_motif_count_dc <- target_motif_count_dc %>% filter(!is.na(Rel_Pos2))



p_dc <- ggplot(target_motif_count_dc,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "reference base", y = "#SNV", fill = "alternative base") +
  theme_minimal() +
  ylim(c(0, 750)) +
  ggtitle("splicing donor creation") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  guides(fill = FALSE)



##########
target_motif_count_ac <- target_motif_count %>% filter(Motif_Type == "splicing acceptor creation")

target_motif_count_ac$Rel_Pos2 <- rep(NA, nrow(target_motif_count_ac))
target_motif_count_ac$Rel_Pos2[target_motif_count_ac$Rel_Pos <= 8] <- 
  9 - target_motif_count_ac$Rel_Pos[target_motif_count_ac$Rel_Pos <= 8] 

target_motif_count_ac$Rel_Pos2[target_motif_count_ac$Rel_Pos >= 9] <- -1


target_motif_count_ac$Rel_Pos2 <- 
  factor(target_motif_count_ac$Rel_Pos2, levels = c("8", "7", "6", "5", "4", "3", "2", "1", "-1"))

target_motif_count_ac <- target_motif_count_ac %>% filter(!is.na(Rel_Pos2))


p_ac <- ggplot(target_motif_count_ac,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "reference base", y = "#SNV", fill = "alternative base") +
  theme_minimal() +
  ylim(c(0, 750)) +
  ggtitle("splicing acceptor creation") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  guides(fill = FALSE)



# p_dc_ac <- plot_grid(p_dc, p_ac, ncol = 2, align = "h", rel_widths = c(0.9, 1))

# plot_grid(p_dc_ac, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.1))


# ggsave("donor_acceptor_creation.png", width = 10, height = 3.5)



p_dd_ad_ac_dc <- plot_grid(p_dd, p_ad, p_dc, p_ac, ncol = 2, align = "h", rel_widths = c(0.9, 1))

plot_grid(p_dd_ad_ac_dc, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.05))

ggsave("snv_motif_dist.png", width = 10, height = 6)



##########
# consistency

edit_dist <- rep(0, nrow(target_motif_count))

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 1 &
            target_motif_count$Ref_Base == "A"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 1 &
            target_motif_count$Alt_Base == "A"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 2 &
            target_motif_count$Ref_Base == "G"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 2 &
            target_motif_count$Alt_Base == "G"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 3 &
            target_motif_count$Ref_Base == "G"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 3 &
            target_motif_count$Alt_Base == "G"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 4 &
            target_motif_count$Ref_Base == "T"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 4 &
            target_motif_count$Alt_Base == "T"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 5 &
            target_motif_count$Ref_Base %in% c("A", "G") &
            target_motif_count$Alt_Base %in% c("C", "T")] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 5 &
            target_motif_count$Alt_Base %in% c("A", "G") & 
            target_motif_count$Ref_Base %in% c("A", "G")] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 6 &
            target_motif_count$Ref_Base == "A"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 6 &
            target_motif_count$Alt_Base == "A"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 7 &
            target_motif_count$Ref_Base == "G"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 7 &
            target_motif_count$Alt_Base == "G"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 8 &
            target_motif_count$Ref_Base == "T"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing donor disruption", "splicing donor creation") &
            target_motif_count$Rel_Pos == 8 &
            target_motif_count$Alt_Base == "T"] <- 1



edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos %in% c(1, 2, 3, 4) &
            target_motif_count$Ref_Base %in% c("C", "T") &
            target_motif_count$Alt_Base %in% c("A", "G")] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos %in% c(1, 2, 3, 4) &
            target_motif_count$Alt_Base %in% c("C", "T") &
            target_motif_count$Ref_Base %in% c("A", "G")] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos == 6 &
            target_motif_count$Ref_Base == "C"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos == 6 &
            target_motif_count$Alt_Base == "C"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos == 7 &
            target_motif_count$Ref_Base == "A"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos == 7 &
            target_motif_count$Alt_Base == "A"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos == 8 &
            target_motif_count$Ref_Base == "G"] <- -1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos %in% c(8, 9) &
            target_motif_count$Alt_Base == "G"] <- 1

edit_dist[target_motif_count$Motif_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") &
            target_motif_count$Rel_Pos %in% c(8, 9) &
            target_motif_count$Ref_Base == "G"] <- -1


target_motif_count$edit_dist <- factor(edit_dist, levels = c(-1, 0, 1), labels = c("get away", "maintain the same dist.", "get close"))

edit_dist_count <- target_motif_count %>% 
  group_by(Motif_Type, edit_dist) %>%
  summarize(ecount = sum(count))


ggplot(edit_dist_count, aes(x = Motif_Type, y = ecount, fill = edit_dist)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  labs(x = "mutation type", y = "#mutation", fill = "edit distance diff.") +
  scale_fill_brewer(palette = "Set3")

ggsave("motif2edit_dist_diff.png", width = 10, height = 3)







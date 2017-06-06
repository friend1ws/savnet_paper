library(dplyr)
library(ggplot2)
library(cowplot)
library(Cairo)

Cairo()

source("../../../conf/plot_config.R")


mes_df <- read.table("../temporary/TCGA.savnet.allele_count.summary.mes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
hb_df <- read.table("../temporary/TCGA.savnet.allele_count.summary.hb.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")


mes_df$splice_class <- factor(mes_df$splice_class,
                              levels = rev(c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                                             "Intron retention", "Complex", "Normal splicing")))

hb_df$splice_class <- factor(hb_df$splice_class,
                             levels = rev(c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                                            "Intron retention", "Complex", "Normal splicing")))


g_mes_d <- ggplot(mes_df %>% filter(splice_class != "Alternative 3'SS" & motif_type == "Donor"), aes(x = splice_class, y = mes_diff, fill = splice_class)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  coord_flip() +
  ggtitle("Donor disruption") +
  ylim(c(-15, 10)) +
  my_theme() +
  theme(axis.title.y = element_text(size = 0)) +
  scale_fill_manual(values = splicing_class_colour) + 
  labs(x = "", y = "") +
  guides(fill = FALSE)

g_mes_a <- ggplot(mes_df %>% filter(splice_class != "Alternative 5'SS" & motif_type == "Acceptor"), aes(x = splice_class, y = mes_diff, fill = splice_class)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.3) +
  coord_flip() +
  ggtitle("Acceptor disruption") +
  labs(x = "", y = "") +
  ylim(c(-15, 10)) +
  my_theme() +
  theme(axis.title.y = element_text(size = 0)) +
  scale_fill_manual(values = splicing_class_colour) + 
  labs(x = "", y = "") +
  guides(fill = FALSE)

ylabel <- ggdraw() + draw_label("Diff. of MaxEnt score", size = 7)

plot_grid(plot_grid(g_mes_d, g_mes_a, ncol = 2, align = "h"), ylabel, ncol = 1, rel_heights = c(1, 0.08))

ggsave("../figure/diff_mes_spliceclass.tiff", width = 16, height = 4.5, dpi = 600, units = "cm")


g_mes_d <- ggplot(mes_df %>% filter(splice_class != "Alternative 3'SS" & motif_type == "Donor"), aes(x = splice_class, y = mes_wt, fill = splice_class)) +
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  coord_flip() +
  ggtitle("Donor disruption") +
  ylim(c(0, 13)) +
  my_theme() +
  theme(plot.margin = unit(c(0.10, 0.10, 0.10, 0.10), "lines"),
        axis.title.x = element_text(size = 0)) +
  scale_fill_manual(values = splicing_class_colour) +
  labs(x = "", y = "") +
  guides(fill = FALSE)

g_mes_a <- ggplot(mes_df %>% filter(splice_class != "Alternative 5'SS" & motif_type == "Acceptor"), aes(x = splice_class, y = mes_wt, fill = splice_class)) +
  geom_boxplot(size = 0.3, outlier.size = 0.3) +
  coord_flip() +
  ggtitle("Acceptor disruption") +
  labs(x = "", y = "") +
  ylim(c(0, 16)) +
  my_theme() +
  theme(plot.margin = unit(c(0.10, 0.10, 0.10, 0.10), "lines"),
        axis.title.x = element_text(size = 0)) +
  scale_fill_manual(values = splicing_class_colour) +
  labs(x = "", y = "") +
  guides(fill = FALSE)

ylabel <- ggdraw() + draw_label("MaxEnt score (before substitution)", size = 7) + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines"))

plot_grid(plot_grid(g_mes_d, g_mes_a, ncol = 2, align = "h"), ylabel, ncol = 1, rel_heights = c(1, 0.12))

ggsave("../figure/mes_wt_spliceclass.tiff", width = 13, height = 4.5, dpi = 600, units = "cm")


pos_colour <- rep("grey30", 9)
pos_colour[4:5] <- "red"

ggplot(mes_df %>% filter(motif_type == "Donor"), 
       aes(x = factor(mut_pos, levels = 1:9), y = mes_diff, fill = is_gsm)) +
  geom_boxplot(size = 0.3, outlier.size = 0.3) +
  labs(fill = "") +
  my_theme() +
  ylim(c(-15, 10)) + 
  labs(x = "Position", y = "Diff. of MaxEnt score", fill = "") +
  scale_x_discrete(limits = 1:9, 
                   labels = c(add_emdash("3"), add_emdash("2"), add_emdash("1"), "+1", "+2", "+3", "+4", "+5", "+6")) +
  scale_fill_manual(values = c("#ef8a62", "#999999")) +
  theme(# axis.text.x = element_text(colour = pos_colour),
        legend.position = "bottom")

ggsave("../figure/diff_mes_mutpos_donor.tiff", width = 7, height = 5, dpi = 600, units = "cm")



pos_colour <- rep("grey30", 7)
pos_colour[5:6] <- "red"

ggplot(mes_df %>% filter(motif_type == "Acceptor"), 
       aes(x = factor(mut_pos, levels = 1:7), y = mes_diff, fill = is_gsm)) +
  geom_boxplot(size = 0.3, outlier.size = 0.3) +
  my_theme() +
  labs(fill = "") +
  ylim(c(-15, 10)) + 
  labs(x = "Position", y = "Diff. of MaxEnt score", fill = "") +
  scale_x_discrete(limits = 1:7, 
                   labels = c("+6", "+5", "+4", "+3", "+2", "+1", add_emdash("1"))) +
  scale_fill_manual(values = c("#ef8a62", "#999999")) +
  theme(# axis.text.x = element_text(colour = pos_colour),
        legend.position = "bottom") 

ggsave("../figure/diff_mes_mutpos_acceptor.tiff", width = 5.5, height = 5, dpi = 600, units = "cm")



ggplot(hb_df %>% filter(splice_class != "Alternative 3'SS"), aes(x = splice_class, y = hb_diff, fill = splice_class)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  coord_flip() +
  ylim(c(-20, 10)) +
  my_theme() +
  scale_fill_manual(values = splicing_class_colour) + 
  labs(x = "", y = "Diff. of H-bond score") +
  guides(fill = FALSE)

ggsave("../figure/diff_hb_spliceclass.tiff", width = 9, height = 4, dpi = 600, units = "cm")


g_mes_d <- ggplot(hb_df %>% filter(splice_class != "Alternative 3'SS"), aes(x = splice_class, y = hb_wt, fill = splice_class)) +
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  coord_flip() +
  ylim(c(5, 25)) +
  ggtitle("Donor disruption") +
  my_theme() +
  theme(plot.margin = unit(c(0.10, 0.10, 0.10, 0.10), "lines"),
        axis.title.x = element_text(size = 0)) +
  scale_fill_manual(values = splicing_class_colour) +
  labs(x = "", y = "") + 
  guides(fill = FALSE)

ylabel <- ggdraw() + draw_label("H-bond score (before substitution)", size = 7) + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines"))

plot_grid(g_mes_d, ylabel, ncol = 1, rel_heights = c(1, 0.12))

ggsave("../figure/hb_wt_spliceclass.tiff", width = 6.5, height = 4.5, dpi = 600, units = "cm")




pos_colour <- rep("grey30", 10)
pos_colour[4:5] <- "red"

ggplot(hb_df, 
       aes(x = factor(mut_pos, levels = 1:9), y = hb_diff, fill = is_gsm)) +
  geom_boxplot(size = 0.3, outlier.size = 0.3) +
  my_theme() +
  labs(fill = "") +
  ylim(c(-20, 10)) + 
  labs(x = "Position", y = "Diff. of H-bond score", fill = "") +
  scale_x_discrete(limits = 1:9, 
                   labels = c(add_emdash("3"), add_emdash("2"), add_emdash("1"), "+1", "+2", "+3", "+4", "+5", "+6")) +
  theme(# axis.text.x = element_text(colour = pos_colour),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#ef8a62", "#999999")) 

ggsave("../figure/diff_hb_mutpos_donor.tiff", width = 7, height = 5, dpi = 600, units = "cm")



library(dplyr)
library(ggplot2)
library(cowplot)

source("../../../conf/plot_config.R")

a <- read.table("../temporary/TCGA.savnet.alt_junc.txt", sep = "\t", header = TRUE, quote = "")

Pos_Diff <- rep(NA, nrow(a))

ind <- a$Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") & a$Exon_Strand == "+"
Pos_Diff[ind] <- a$Exon_Start[ind] - a$Junc_Pos[ind] 

ind <- a$Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") & a$Exon_Strand == "-"
Pos_Diff[ind] <- a$Junc_Pos[ind] - a$Exon_End[ind] 

ind <- a$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation") & a$Exon_Strand == "+"
Pos_Diff[ind] <- a$Junc_Pos[ind] - a$Exon_End[ind] 

ind <- a$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation") & a$Exon_Strand == "-"
Pos_Diff[ind] <- a$Exon_Start[ind] - a$Junc_Pos[ind] 


a$Pos_Diff <- Pos_Diff


a <- a %>% filter(Pos_Diff <= 100) %>% filter(Pos_Diff >= -300)
a_donor <- a %>% filter(Mutation_Type == "splicing donor creation")
a_acceptor <- a %>% filter(Mutation_Type == "splicing acceptor creation")

a$Mutation_Type <- factor(a$Mutation_Type, 
                          levels = c("splicing donor disruption", "splicing acceptor disruption", 
                            "splicing donor creation", "splicing acceptor creation"))

a$Splicing_Class <- rep("Alternative 5'SS", nrow(a))
a$Splicing_Class[a$Mutation_Type %in% c("splicing acceptor creation", "splicing acceptor disruption")] <- "Alternative 3'SS"



p_dd <- ggplot(a %>% filter(Mutation_Type == "splicing donor disruption"),
               aes(x = Pos_Diff, fill = Splicing_Class)) + 
  geom_histogram(binwidth = 5, colour = "grey30", size = 0.2) +
  my_theme() +
  theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "lines")) + 
  xlim(c(-300, 100)) +
  ggtitle("Donor disruption") +
  labs(x = "", y = "") +
  scale_fill_manual(values = splicing_class_colour) + 
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(size = 0.3, xintercept = 0, colour="#d73027", linetype = "longdash") +
  guides(fill = FALSE)

p_dc <- ggplot(a %>% filter(Mutation_Type == "splicing donor creation"),
               aes(x = Pos_Diff, fill = Splicing_Class)) + 
  geom_histogram(binwidth = 5, colour = "grey30", size = 0.2) +
  my_theme() +
  theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "lines")) + 
  xlim(c(-300, 100)) +
  ggtitle("Donor creation") +
  labs(x = "", y = "") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(size = 0.3, xintercept = 0, colour="#d73027", linetype = "longdash") +
  guides(fill = FALSE)

p_ad <- ggplot(a %>% filter(Mutation_Type == "splicing acceptor disruption"),
               aes(x = Pos_Diff, fill = Splicing_Class)) + 
  geom_histogram(binwidth = 5, colour = "grey30", size = 0.2) +
  my_theme() +
  theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "lines")) +
  scale_x_reverse(limits=c(100, -300)) +
  ggtitle("Acceptor disruption") +
  labs(x = "", y = "") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(size = 0.3, xintercept = 0, colour="#d73027", linetype = "longdash") +
  guides(fill = FALSE)

p_ac <- ggplot(a %>% filter(Mutation_Type == "splicing acceptor creation"),
               aes(x = Pos_Diff, fill = Splicing_Class)) + 
  geom_histogram(binwidth = 5, colour = "grey30", size = 0.2) +
  my_theme() +
  theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "lines")) +
  scale_x_reverse(limits=c(100, -300)) +
  ggtitle("Acceptor creation") +
  labs(x = "", y = "") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(size = 0.3, xintercept = 0, colour="#d73027", linetype = "longdash") +
  guides(fill = FALSE) 

# p_dc_ac <- plot_grid(p_dc, p_ac, ncol = 2, align = "h")
# p_dd_ad <- plot_grid(p_dd, p_ad, ncol = 2, align = "h")

# p_dc_ac_dd_ad <- plot_grid(p_dc_ac, p_dd_ad, ncol = 1, align = "v")

p_dc_ac_dd_ad <- plot_grid(p_dc, p_ac, p_dd, p_ad, ncol = 2, align = "hv")

xlabel <- ggdraw() + draw_label("Position", size = 7)
ylabel <- ggdraw() + draw_label("Abnormal splicing event count", angle = 90, size = 7)

p_dc_ac_dd_ad_xl <- plot_grid(p_dc_ac_dd_ad, xlabel, ncol = 1, align = "v", rel_heights = c(2, 0.1))

plot_grid(ylabel, p_dc_ac_dd_ad_xl, ncol = 2, align = "h", rel_widths = c(0.04, 1))


ggsave("../figure/alt_junc_pos.tiff", width = 10, height = 7, dpi = 600, units = "cm")



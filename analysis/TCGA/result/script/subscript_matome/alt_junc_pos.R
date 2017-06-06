library(dplyr)
library(ggplot2)
library(cowplot)
library(Cairo)

Cairo()

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
  theme(plot.margin = unit(c(0.21, 0.21, 0.21, 0.21), "lines"),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) + 
  # xlim(c(-300, 100)) +
  # scale_x_continuous(limits = c(-300, 100), breaks = c(-300, -200, -100, 0, 100), labels = c("-300", "-200", "-100", "-1 +1", "+100")) +
  scale_x_continuous(limits = c(-200, 100), breaks = c(-200, -100, 0, 100), 
                     labels = c(add_emdash("200"), add_emdash("100"), paste(add_emdash("1"), "+1", sep = " "), "+100")) + 
  ggtitle("Cryptic 5'SS by \ndonor disruption") +
  labs(x = "", y = "") +
  scale_fill_manual(values = splicing_class_colour) + 
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(size = 0.3, xintercept = 0, colour="#d73027", linetype = "longdash") +
  guides(fill = FALSE)

p_dc <- ggplot(a %>% filter(Mutation_Type == "splicing donor creation"),
               aes(x = Pos_Diff, fill = Splicing_Class)) + 
  geom_histogram(binwidth = 5, colour = "grey30", size = 0.2) +
  my_theme() +
  theme(plot.margin = unit(c(0.21, 0.21, 0.21, 0.21), "lines"),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) + 
  # xlim(c(-300, 100)) +
  # scale_x_continuous(limits = c(-300, 100), breaks = c(-300, -200, -100, 0, 100), labels = c("-300", "-200", "-100", "-1 +1", "+100")) + 
  scale_x_continuous(limits = c(-200, 100), breaks = c(-200, -100, 0, 100), 
                     labels = c(add_emdash("200"), add_emdash("100"), paste(add_emdash("1"), "+1", sep = " "), "+100")) + 
  ggtitle("Alternative 5'SS by \ndonor creation") +
  labs(x = "", y = "") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(size = 0.3, xintercept = 0, colour="#d73027", linetype = "longdash") +
  guides(fill = FALSE)

p_ad <- ggplot() +
  geom_rect(aes(xmin=25, xmax=5, ymin=0, ymax=Inf), fill = "#fccde5", alpha = 0.70) +
  geom_histogram(data = a %>% filter(Mutation_Type == "splicing acceptor disruption"),
                 aes(x = Pos_Diff, fill = Splicing_Class), binwidth = 5, colour = "grey30", size = 0.2) +
  my_theme() +
  theme(plot.margin = unit(c(0.21, 0.21, 0.21, 0.21), "lines"),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0)) +
  # scale_x_reverse(limits=c(100, -300)) +
  # scale_x_continuous(trans = "reverse", limits = c(100, -300), breaks = c(100, 0, -100, -200, -300), labels = c("+100", "+1 -1", "-100", "-200", "-300")) +
  scale_x_continuous(trans = "reverse", limits = c(100, -200), breaks = c(100, 0, -100, -200), 
                     labels = c("+100", paste("+1", add_emdash("1"), sep = " "), add_emdash("100"), add_emdash("200"))) +
  ggtitle("Cryptic 3'SS by \nacceptor disruption") +
  labs(x = "", y = "") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(size = 0.3, xintercept = 0, colour="#d73027", linetype = "longdash") +
  guides(fill = FALSE)

p_ac <- ggplot() +
  geom_rect(aes(xmin=25, xmax=5, ymin=0, ymax=Inf), fill = "#fccde5", alpha = 0.70) + 
  geom_histogram(data = a %>% filter(Mutation_Type == "splicing acceptor creation"),
                 aes(x = Pos_Diff, fill = Splicing_Class), binwidth = 5, colour = "grey30", size = 0.2) +
  my_theme() +
  theme(plot.margin = unit(c(0.21, 0.21, 0.21, 0.21), "lines"),
        axis.title.x = element_text(size = 0), 
        axis.title.y = element_text(size = 0)) +
  # scale_x_reverse(limits=c(100, -300)) +
  # scale_x_continuous(trans = "reverse", limits = c(100, -300), breaks = c(100, 0, -100, -200, -300), labels = c("+100", "+1 -1", "-100", "-200", "-300")) +
  scale_x_continuous(trans = "reverse", limits = c(100, -200), breaks = c(100, 0, -100, -200), 
                     labels = c("+100", paste("+1", add_emdash("1"), sep = " "), add_emdash("100"), add_emdash("200"))) +
  ggtitle("Alternative 3'SS by \nacceptor creation") +
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

p_dc_ac_dd_ad_xl <- plot_grid(p_dc_ac_dd_ad, ggdraw() + draw_label("", size = 7), xlabel, ncol = 1, align = "v", rel_heights = c(2, 0.15, 0.1))
# p_dc_ac_dd_ad_xl <- plot_grid(p_dc_ac_dd_ad, xlabel, ncol = 1, align = "v", rel_heights = c(2, 0.1))

plot_grid(ylabel, p_dc_ac_dd_ad_xl, ncol = 2, align = "h", rel_widths = c(0.05, 1))


ggsave("../figure/alt_junc_pos.tiff", width = 10, height = 7, dpi = 600, units = "cm")



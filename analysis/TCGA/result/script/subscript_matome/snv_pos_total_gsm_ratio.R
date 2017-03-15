library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("../../../conf/plot_config.R")

gg_color_hue6 <- hcl(h = seq(15, 375, length = 7), l=65, c=100)[1:6]
base_col <- c("A" = gg_color_hue6[3], "C" = gg_color_hue6[5], "G" = gg_color_hue6[2], "T" = gg_color_hue6[1])

get_pos_df <- function(motif_pos) {
  
  motif_pos_sp1 <- strsplit(motif_pos, ":")
  tchr <- paste("chr", unlist(lapply(motif_pos_sp1, '[', 1)), sep = "")
  
  motif_pos_sp2 <- strsplit(unlist(lapply(motif_pos_sp1, '[', 2)), ",")
  tstrand <- unlist(lapply(motif_pos_sp2, '[', 2))
  
  motif_pos_sp3 <- strsplit(unlist(lapply(motif_pos_sp2, '[', 1)), "-")
  tstart <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 1)))
  tend <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 2)))
  
  return(data.frame(chr = tchr, start = tstart, end = tend, strand = tstrand))
  
}

snv_info <- read.table("../temporary/TCGA.savnet.allele_count.summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "") %>%
  filter(Ref_Mut != "-" & Alt_Mut != "-")

snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "A"] <- "S"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "T"] <- "A"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "S"] <- "T"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "C"] <- "S"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "G"] <- "C"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "S"] <- "G"

snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "A"] <- "S"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "T"] <- "A"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "S"] <- "T"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "C"] <- "S"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "G"] <- "C"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "S"] <- "G"


snv_info$Is_GSM <- "no"
snv_info$Is_GSM[snv_info$GenomonSplicingMutation != "---"] <- "yes"




snv_motif_count <- snv_info %>% 
  group_by(Rel_Start_Motif, Ref_Mut, Alt_Mut, Type_Motif, Is_GSM) %>% 
  summarize(count = n())

snv_motif_count$Rel_Pos <- snv_motif_count$Rel_Start_Motif
snv_motif_count$Ref_Base <- snv_motif_count$Ref_Mut
snv_motif_count$Alt_Base <- snv_motif_count$Alt_Mut
snv_motif_count$Mutation_Type <- rep("splicing donor disruption", nrow(snv_motif_count))
snv_motif_count$Mutation_Type[snv_motif_count$Type_Motif == "acceptor"] <- "splicing acceptor disruption"







##########
exon_size_d <- 3

snv_motif_count_dd <- snv_motif_count %>% filter(Type_Motif == "donor")

snv_motif_count_dd$Rel_Start_Motif2 <- rep(NA, nrow(snv_motif_count_dd))
snv_motif_count_dd$Rel_Start_Motif2[snv_motif_count_dd$Rel_Start_Motif <= exon_size_d] <- 
  snv_motif_count_dd$Rel_Start_Motif[snv_motif_count_dd$Rel_Start_Motif <= exon_size_d] - exon_size_d - 1

snv_motif_count_dd$Rel_Start_Motif2[snv_motif_count_dd$Rel_Start_Motif > exon_size_d] <- 
  snv_motif_count_dd$Rel_Start_Motif[snv_motif_count_dd$Rel_Start_Motif > exon_size_d] - exon_size_d


p_dd_total <- ggplot(snv_motif_count_dd,
               aes(x = Ref_Mut, y = count, fill = Alt_Mut)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif2, nrow = 1) +
  labs(x = "", y = "Total mutation count", fill = "Alternative base") +
  my_theme() +
  ggtitle("Donor Disruption") +
  theme(# axis.text.x = element_text(size = rel(1)),
        # axis.text.y = element_text(size = rel(1)),
        # axis.title = element_text(size = rel(1)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  scale_y_continuous(limits = c(0, 16000)) + 
  guides(fill = FALSE)


snv_gsm_ratio_dd <- snv_motif_count_dd %>% 
  group_by(Rel_Start_Motif2, Ref_Mut, Is_GSM) %>% 
  summarize(sum = sum(count)) %>% 
  filter(sum >= 10) %>%
  spread(key = Is_GSM, value = sum, fill = 0) 

snv_gsm_ratio_dd$ratio_lower <- nrow(snv_gsm_ratio_dd)
snv_gsm_ratio_dd$ratio_median <- nrow(snv_gsm_ratio_dd)
snv_gsm_ratio_dd$ratio_upper <- nrow(snv_gsm_ratio_dd)

for (i in 1:nrow(snv_gsm_ratio_dd)) {
  snv_gsm_ratio_dd$ratio_lower[i] <- qbeta(0.05, snv_gsm_ratio_dd$yes[i] + 1, snv_gsm_ratio_dd$no[i] + 1)
  snv_gsm_ratio_dd$ratio_median[i] <- qbeta(0.5, snv_gsm_ratio_dd$yes[i] + 1, snv_gsm_ratio_dd$no[i] + 1)
  snv_gsm_ratio_dd$ratio_upper[i] <- min(qbeta(0.95, snv_gsm_ratio_dd$yes[i] + 1, snv_gsm_ratio_dd$no[i] + 1), 0.25)
}


p_dd_gsm <- ggplot() +
  geom_bar(data = snv_motif_count_dd %>% filter(Is_GSM == "yes"),
           aes(x = Ref_Mut, y = count, fill = Alt_Mut), stat = "identity") +
  geom_point(data = snv_gsm_ratio_dd, aes(x = Ref_Mut, y = 12000 * ratio_median), colour = "#7570b3", size = 0.4, alpha = 0.8) +
  geom_segment(data = snv_gsm_ratio_dd, aes(x = Ref_Mut, xend = Ref_Mut,
                                            y = 12000 * ratio_lower, yend = 12000 * ratio_upper), colour = "#7570b3", alpha = 0.4) +
  # geom_line(data = snv_gsm_ratio, aes(x = Ref_Mut, y = 12000 * ratio), colour = "#7570b3", alpha = 0.8) +
  facet_wrap( ~ Rel_Start_Motif2, nrow = 1) +
  labs(x = "Reference base", y = "SAV count", fill = "Alternative base") +
  my_theme() +
  theme(# axis.text.x = element_text(size = rel(1)),
        # axis.text.y = element_text(size = rel(1)),
        # axis.title = element_text(size = rel(1)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  scale_y_continuous(limits = c(0, 3000), sec.axis = sec_axis(~ . * (1 / 12000), name = "")) +
  guides(fill = FALSE)



##########
# splicing acceptor motif
intron_size_a <- 6

snv_motif_count_ad <- snv_motif_count %>% filter(Type_Motif == "acceptor")

snv_motif_count_ad$Rel_Start_Motif2 <- rep(NA, nrow(snv_motif_count_ad))
snv_motif_count_ad$Rel_Start_Motif2[snv_motif_count_ad$Rel_Start_Motif <= intron_size_a] <- 
  intron_size_a - snv_motif_count_ad$Rel_Start_Motif[snv_motif_count_ad$Rel_Start_Motif <= intron_size_a] + 1

snv_motif_count_ad$Rel_Start_Motif2[snv_motif_count_ad$Rel_Start_Motif > intron_size_a] <- 
  intron_size_a - snv_motif_count_ad$Rel_Start_Motif[snv_motif_count_ad$Rel_Start_Motif > intron_size_a]  


snv_motif_count_ad$Rel_Start_Motif2 <- factor(snv_motif_count_ad$Rel_Start_Motif2,
                                              levels = rev(c(-1, 1, 2, 3, 4, 5, 6)),
                                              labels = rev(c("-1", "1", "2", "3", "4", "5", "6")))


p_ad_total <- ggplot(snv_motif_count_ad,
                     aes(x = Ref_Mut, y = count, fill = Alt_Mut)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif2, nrow = 1) +
  labs(x = "", y = "", fill = "Alternative base") +
  my_theme() +
  ggtitle("Acceptor disruption") +
  theme(# axis.text.x = element_text(size = rel(1)),
        # axis.text.y = element_text(size = rel(1)),
        # axis.title = element_text(size = rel(1)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  scale_y_continuous(limits = c(0, 16000)) + 
  guides(fill = FALSE)


snv_gsm_ratio_ad <- snv_motif_count_ad %>% 
  group_by(Rel_Start_Motif2, Ref_Mut, Is_GSM) %>% 
  summarize(sum = sum(count)) %>% 
  filter(sum >= 10) %>%
  spread(key = Is_GSM, value = sum, fill = 0) 

snv_gsm_ratio_ad$ratio_lower <- nrow(snv_gsm_ratio_ad)
snv_gsm_ratio_ad$ratio_median <- nrow(snv_gsm_ratio_ad)
snv_gsm_ratio_ad$ratio_upper <- nrow(snv_gsm_ratio_ad)

for (i in 1:nrow(snv_gsm_ratio_ad)) {
  snv_gsm_ratio_ad$ratio_lower[i] <- qbeta(0.05, snv_gsm_ratio_ad$yes[i] + 1, snv_gsm_ratio_ad$no[i] + 1)
  snv_gsm_ratio_ad$ratio_median[i] <- qbeta(0.5, snv_gsm_ratio_ad$yes[i] + 1, snv_gsm_ratio_ad$no[i] + 1)
  snv_gsm_ratio_ad$ratio_upper[i] <- min(qbeta(0.95, snv_gsm_ratio_ad$yes[i] + 1, snv_gsm_ratio_ad$no[i] + 1), 0.25)
}

p_ad_gsm <- ggplot() +
  geom_bar(data = snv_motif_count_ad %>% filter(Is_GSM == "yes"),
                   aes(x = Ref_Mut, y = count, fill = Alt_Mut), stat = "identity") +
  geom_point(data = snv_gsm_ratio_ad, aes(x = Ref_Mut, y = 12000 * ratio_median), colour = "#7570b3", size = 0.4, alpha = 0.8) +
  geom_segment(data = snv_gsm_ratio_ad, aes(x = Ref_Mut, xend = Ref_Mut,
                                            y = 12000 * ratio_lower, yend = 12000 * ratio_upper), colour = "#7570b3", alpha = 0.4) +
  facet_wrap( ~ Rel_Start_Motif2, nrow = 1) +
  labs(x = "Reference base", y = "", fill = "Alternative base") +
  my_theme() +
  theme(# axis.text.x = element_text(size = rel(1)),
        # axis.text.y = element_text(size = rel(1)),
        # axis.title = element_text(size = rel(1)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  scale_y_continuous(limits = c(0, 3000), sec.axis = sec_axis(~ . * (1 / 12000), name = "Splicing ratio")) +
  guides(fill = FALSE)


##########
# legend

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p_dummy_for_legend <- ggplot(snv_motif_count %>% filter(Type_Motif == "donor"),
                             aes(x = Ref_Mut, y = count, fill = Alt_Mut)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif, nrow = 1) +
  labs(x = "Reference base", y = "#SNV", fill = "Alternative base") +
  my_theme() +
  ggtitle("Donor") +
  theme(# axis.text.x = element_text(size = rel(1)),
        # axis.text.y = element_text(size = rel(1)),
        # axis.title = element_text(size = rel(1)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) 


# p_total <- plot_grid(p_dd_total, p_ad_total, ncol = 2, align = "h", rel_widths = c(1, 0.9))
# p_gsm <- plot_grid(p_dd_gsm, p_ad_gsm, ncol = 2, align = "h", rel_widths = c(1, 0.9))
# p_ratio <- plot_grid(p_dd_ratio, p_ad_ratio, ncol = 2, align = "h", rel_widths = c(1, 0.9))

p_donor <- plot_grid(p_dd_total, p_dd_gsm, ncol = 1, rel_heights = c(1.1, 1), align = "v")
p_acceptor <- plot_grid(p_ad_total, p_ad_gsm, ncol = 1, rel_heights = c(1.1, 1), align = "v")

plot_grid(plot_grid(p_donor, p_acceptor, ncol = 2, align = "h", rel_widths = c(1, 0.85)), g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(2.2, 0.1))


ggsave("../figure/snv_pos_total_gsm_ratio.tiff", width = 20, height = 9, dpi = 600, units = "cm")




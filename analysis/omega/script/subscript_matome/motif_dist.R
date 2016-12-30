library(dplyr)
library(ggplot2)
library(cowplot)

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

snv_info <- read.table("../matome/omega.motif_summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# snv_info <- read.table("omega.motif_summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

snv_motif_count <- snv_info %>% 
  group_by(Rel_Pos, Ref_Base, Alt_Base, Mutation_Type) %>% 
  summarize(count = n())


##########
# get exon intron sizes

canonical_count_dd <- snv_info  %>% 
  filter(Mutation_Type == "splicing donor disruption", Is_Canonical == "canonical") %>% 
  group_by(Rel_Pos) %>% summarize(count = n())

if (length(canonical_count_dd$Rel_Pos) != 2) 
  stop("canonical info of donor splicing motif is inconsistent!")

snv_info_d <- snv_info %>% 
  filter(Mutation_Type %in% c("splicing donor disruption", "splicing donor creation"))

pos_df_d <- get_pos_df(snv_info_d$Motif_Pos)

motif_len_d <- unique(pos_df_d$end - pos_df_d$start + 1)
if (length(motif_len_d) != 1) stop("donor motif size is inconsistent!")

exon_size_d <- min(canonical_count_dd$Rel_Pos) - 1
intron_size_d <- motif_len_d - exon_size_d


motif_len_d <- unique(pos_df_d$end - pos_df_d$start + 1)
if (length(motif_len_d) != 1) stop("donor motif size is inconsistent!")

canonical_count_ad <- snv_info  %>% 
  filter(Mutation_Type == "splicing acceptor disruption", Is_Canonical == "canonical") %>% 
  group_by(Rel_Pos) %>% summarize(count = n())

if (length(canonical_count_ad$Rel_Pos) != 2) 
  stop("canonical info of acceptor splicing motif is inconsistent!")


snv_info_a <- snv_info %>% 
  filter(Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation"))

pos_df_a <- get_pos_df(snv_info_a$Motif_Pos)

motif_len_a <- unique(pos_df_a$end - pos_df_a$start + 1)
if (length(motif_len_a) != 1) stop("donor motif size is inconsistent!")


intron_size_a <- max(canonical_count_ad$Rel_Pos) 
exon_size_a <- motif_len_a - intron_size_a



##########
# splicing donor motif (target)
# pos_colour <- rep("grey30", 8)
# pos_colour[3:4] <- "red"

snv_motif_count_dd <- snv_motif_count %>% filter(Mutation_Type == "splicing donor disruption")

snv_motif_count_dd$Rel_Pos2 <- rep(NA, nrow(snv_motif_count_dd))
snv_motif_count_dd$Rel_Pos2[snv_motif_count_dd$Rel_Pos <= exon_size_d] <- 
  snv_motif_count_dd$Rel_Pos[snv_motif_count_dd$Rel_Pos <= exon_size_d] - exon_size_d - 1

snv_motif_count_dd$Rel_Pos2[snv_motif_count_dd$Rel_Pos > exon_size_d] <- 
  snv_motif_count_dd$Rel_Pos[snv_motif_count_dd$Rel_Pos > exon_size_d] - exon_size_d


snv_motif_count_dd$Rel_Pos2 <- 
  factor(snv_motif_count_dd$Rel_Pos2, levels = unique(as.character(sort(snv_motif_count_dd$Rel_Pos2))))

# pos_colour <- rep("grey30", motif_len_d)
# pos_colour[canonical_count_dd$Rel_Pos] <- "red"

p_dd <- ggplot(snv_motif_count_dd,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "Reference base", y = "#SNV", fill = "Alternative base") +
  ylim(c(0, 3000)) +
  theme_minimal() +
  ggtitle("Donor disruption") +
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


snv_motif_count_ad <- snv_motif_count %>% filter(Mutation_Type == "splicing acceptor disruption")

snv_motif_count_ad$Rel_Pos2 <- rep(NA, nrow(snv_motif_count_ad))
snv_motif_count_ad$Rel_Pos2[snv_motif_count_ad$Rel_Pos <= intron_size_a] <- 
  intron_size_a - snv_motif_count_ad$Rel_Pos[snv_motif_count_ad$Rel_Pos <= intron_size_a] + 1

snv_motif_count_ad$Rel_Pos2[snv_motif_count_ad$Rel_Pos > intron_size_a] <- 
  intron_size_a - snv_motif_count_ad$Rel_Pos[snv_motif_count_ad$Rel_Pos > intron_size_a]  


snv_motif_count_ad$Rel_Pos2 <- 
  factor(snv_motif_count_ad$Rel_Pos2, levels = rev(unique(as.character(sort(snv_motif_count_ad$Rel_Pos2)))))


p_ad <- ggplot(snv_motif_count_ad,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Rel_Pos2, nrow = 1) +
  labs(x = "Reference base", y = "#SNV", fill = "Alternative base") +
  ylim(c(0, 3000)) +
  theme_minimal() +
  ggtitle("Acceptor disruption") +
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

p_dummy_for_legend <- ggplot(snv_motif_count_dd,
                             aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "Reference base", y = "#SNV", fill = "Alternative base") +
  theme_minimal() +
  ggtitle("Donor disruption") +
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
snv_motif_count_dc <- snv_motif_count %>% filter(Mutation_Type == "splicing donor creation")

snv_motif_count_dc$Rel_Pos2 <- rep(NA, nrow(snv_motif_count_dc))
snv_motif_count_dc$Rel_Pos2[snv_motif_count_dc$Rel_Pos <= exon_size_d] <- 
  snv_motif_count_dc$Rel_Pos[snv_motif_count_dc$Rel_Pos <= exon_size_d] - exon_size_d - 1

snv_motif_count_dc$Rel_Pos2[snv_motif_count_dc$Rel_Pos > exon_size_d] <- 
  snv_motif_count_dc$Rel_Pos[snv_motif_count_dc$Rel_Pos > exon_size_d] - exon_size_d


snv_motif_count_dc$Rel_Pos2 <- 
  factor(snv_motif_count_dc$Rel_Pos2, levels = unique(as.character(sort(snv_motif_count_dc$Rel_Pos2))))


p_dc <- ggplot(snv_motif_count_dc,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1, drop = FALSE) +
  labs(x = "Reference base", y = "#SNV", fill = "Alternative base") +
  theme_minimal() +
  ylim(c(0, 750)) +
  ggtitle("Donor creation") +
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
snv_motif_count_ac <- snv_motif_count %>% filter(Mutation_Type == "splicing acceptor creation")

snv_motif_count_ac$Rel_Pos2 <- rep(NA, nrow(snv_motif_count_ac))
snv_motif_count_ac$Rel_Pos2[snv_motif_count_ac$Rel_Pos <= intron_size_a] <- 
  intron_size_a - snv_motif_count_ac$Rel_Pos[snv_motif_count_ac$Rel_Pos <= intron_size_a] + 1

snv_motif_count_ac$Rel_Pos2[snv_motif_count_ac$Rel_Pos > intron_size_a] <- 
  intron_size_a - snv_motif_count_ac$Rel_Pos[snv_motif_count_ac$Rel_Pos > intron_size_a]  




snv_motif_count_ac$Rel_Pos2 <- 
  factor(snv_motif_count_ac$Rel_Pos2, levels = rev(unique(as.character(sort(snv_motif_count_ac$Rel_Pos2)))))


p_ac <- ggplot(snv_motif_count_ac,
       aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "Reference base", y = "#SNV", fill = "Alternative base") +
  theme_minimal() +
  ylim(c(0, 750)) +
  ggtitle("Acceptor creation") +
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

ggsave("../matome/snv_motif_dist.png", width = 10, height = 6)



##########
# consistency

if (FALSE) {

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

ggsave("../matome/motif2edit_dist_diff.png", width = 10, height = 3)


write.table(edit_dist_count, "../matome/motif2edit_dist_diff.txt",
             quote = FALSE, row.names = FALSE, sep = "\t")

}





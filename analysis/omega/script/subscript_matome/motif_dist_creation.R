library(dplyr)
library(ggplot2)
library(cowplot)

source("subscript_matome/plot_config.R")

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
  labs(x = "Reference base", y = "SAV count", fill = "Alternative base") +
  my_theme() +
  ggtitle("Donor creation") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 750)) +
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
  labs(x = "Reference base", y = "SAV count", fill = "Alternative base") +
  my_theme() +
  ggtitle("Acceptor creation") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 250)) +
  guides(fill = FALSE)


# legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p_dummy_for_legend <- ggplot(snv_motif_count_dc,
                             aes(x = Ref_Base, y = count, fill = Alt_Base)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Pos2, nrow = 1) +
  labs(x = "Reference base", y = "SAV count", fill = "Alternative base") +
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


p_ac_dc <- plot_grid(p_dc, p_ac, align = "h", rel_widths = c(1, 0.9))

plot_grid(p_ac_dc, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.1))

ggsave("../matome/snv_motif_dist_creation.pdf", width = 10, height = 3.3)




library(dplyr)
library(tidyr)
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

snv_info <- read.table("../output/omega.splicing_mutation.info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "") %>%
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

p_dd_total <- ggplot(snv_motif_count %>% filter(Type_Motif == "donor"),
               aes(x = Ref_Mut, y = count, fill = Alt_Mut)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif, nrow = 1) +
  labs(x = "Reference base", y = "Total mutation count", fill = "Alternative base") +
  theme_minimal() +
  ggtitle("Donor") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  guides(fill = FALSE)

p_dd_gsm <- ggplot(snv_motif_count %>% filter(Type_Motif == "donor", Is_GSM == "yes"),
                     aes(x = Ref_Mut, y = count, fill = Alt_Mut)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif, nrow = 1) +
  labs(x = "Reference base", y = "SASM count", fill = "Alternative base") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  guides(fill = FALSE)


snv_gsm_ratio <- snv_motif_count %>% 
  filter(Type_Motif == "donor") %>% 
  group_by(Rel_Start_Motif, Ref_Mut, Is_GSM) %>% 
  summarize(sum = sum(count)) %>% 
  spread(key = Is_GSM, value = sum, fill = 0) %>% 
  mutate(ratio = yes / (yes + no))


p_dd_ratio <- ggplot(snv_gsm_ratio,
                   aes(x = Ref_Mut, y = ratio)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif, nrow = 1) +
  labs(x = "Reference base", y = "Splicing ratio", fill = "Alternative base") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") 

##########
# splicing acceptor motif

p_ad_total <- ggplot(snv_motif_count %>% filter(Type_Motif == "acceptor"),
                     aes(x = Ref_Mut, y = count, fill = Alt_Mut)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif, nrow = 1) +
  labs(x = "Reference base", y = "Total mutation count", fill = "Alternative base") +
  theme_minimal() +
  ggtitle("Acceptor") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  guides(fill = FALSE)

p_ad_gsm <- ggplot(snv_motif_count %>% filter(Type_Motif == "acceptor", Is_GSM == "yes"),
                   aes(x = Ref_Mut, y = count, fill = Alt_Mut)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif, nrow = 1) +
  labs(x = "Reference base", y = "SASM count", fill = "Alternative base") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) +
  guides(fill = FALSE)


snv_gsm_ratio <- snv_motif_count %>% 
  filter(Type_Motif == "acceptor") %>% 
  group_by(Rel_Start_Motif, Ref_Mut, Is_GSM) %>% 
  summarize(sum = sum(count)) %>% 
  spread(key = Is_GSM, value = sum, fill = 0) %>% 
  mutate(ratio = yes / (yes + no))


p_ad_ratio <- ggplot(snv_gsm_ratio,
                     aes(x = Ref_Mut, y = ratio)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Rel_Start_Motif, nrow = 1) +
  labs(x = "Reference base", y = "Splicing ratio", fill = "Alternative base") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") 


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
  theme_minimal() +
  ggtitle("Donor") +
  theme(axis.text.x = element_text(size = rel(1)),
        axis.text.y = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = base_col) 


p_total <- plot_grid(p_dd_total, p_ad_total, ncol = 2, align = "h", rel_widths = c(1, 0.9))
p_gsm <- plot_grid(p_dd_gsm, p_ad_gsm, ncol = 2, align = "h", rel_widths = c(1, 0.9))
p_ratio <- plot_grid(p_dd_ratio, p_ad_ratio, ncol = 2, align = "h", rel_widths = c(1, 0.9))


plot_grid(p_total, p_gsm, g_legend(p_dummy_for_legend), p_ratio, ncol = 1, rel_heights = c(1.2, 1, 0.1, 1))


ggsave("../output/snv_pos_total_gsm_ratio.pdf", width = 10, height = 10)




library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("../../../conf/plot_config.R")

splicing_mut_info <- read.table("../temporary/TCGA.savnet.allele_count.summary.txt", header = TRUE, sep = "\t", 
                                as.is=TRUE, quote="", stringsAsFactors = FALSE)

                                
splicing_mut_info_filt <- splicing_mut_info 
##########

##########
# organize the splicing type information
splicing_mut_info_filt[
  grep(";", splicing_mut_info_filt$GenomonSplicingMutation),
  "GenomonSplicingMutation"] <- "complex"

splicing_mut_info_filt[
  splicing_mut_info_filt$GenomonSplicingMutation == "---",
  "GenomonSplicingMutation"] <- "no-change"

splicing_mut_info_filt[
  splicing_mut_info_filt$GenomonSplicingMutation == "opposite-side-intron-retention",
  "GenomonSplicingMutation"] <- "intron-retention"

splicing_mut_info_filt[
  splicing_mut_info_filt$GenomonSplicingMutation == "intronic-alternative-5'-splice-site",
  "GenomonSplicingMutation"] <- "alternative-5'-splice-site"

splicing_mut_info_filt[
  splicing_mut_info_filt$GenomonSplicingMutation == "intronic-alternative-3'-splice-site",
  "GenomonSplicingMutation"] <- "alternative-3'-splice-site"

splicing_mut_info_filt$GenomonSplicingMutation <-
  factor(splicing_mut_info_filt$GenomonSplicingMutation,
         levels = c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
                    "intron-retention", "complex", "no-change"),
         labels = c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                    "Intron retention", "Complex", "Normal splicing"))
 


##########

splicing_mut_info_filt_snv <- splicing_mut_info_filt %>% 
  filter(Ref_Mut != "-", Alt_Mut != "-")


##########
# obtain the ratio
splicing_mut_info_filt_snv_count <- splicing_mut_info_filt_snv %>% 
  group_by(Rel_Start_Motif, Type_Motif, GenomonSplicingMutation) %>% 
  summarize(count = n())

splicing_mut_info_filt_snv_count_total <- splicing_mut_info_filt_snv_count %>% 
  group_by(Rel_Start_Motif, Type_Motif) %>% 
  summarize(total_count = sum(count))

splicing_mut_info_filt_snv_ratio <- 
  left_join(splicing_mut_info_filt_snv_count, splicing_mut_info_filt_snv_count_total, by = c("Rel_Start_Motif", "Type_Motif")) %>%
  mutate(ratio = count / total_count) 

write.table(splicing_mut_info_filt_snv_ratio, 
  "../temporary/splicing_snv_ratio.txt", quote = FALSE, row.names = FALSE, sep = "\t")

 
pos_colour <- rep("grey30", 10)
pos_colour[4:5] <- "red"

p_donor_count <- ggplot(splicing_mut_info_filt_snv_count %>% 
                          filter(Type_Motif == "donor", GenomonSplicingMutation != "Normal splicing"), 
                        aes(x = Rel_Start_Motif, y = count, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", y = "SAV count", fill = "") +
  ggtitle("Donor") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(# axis.text.x = element_text(colour = pos_colour),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:9, 
                   # labels = c("M", "A", "G", "G", "T", "R", "A", "G", "T")) +
                   labels = c("-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5", "+6")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 3000)) + 
  guides(fill = FALSE)


p_donor_ratio <- ggplot(splicing_mut_info_filt_snv_ratio %>% 
         filter(Type_Motif == "donor", GenomonSplicingMutation != "Normal splicing"), 
       aes(x = Rel_Start_Motif, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", y = "SAV ratio", fill = "") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(# axis.text.x = element_text(colour = pos_colour),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:9, 
                   # labels = c("M", "A", "G", "G", "T", "R", "A", "G", "T")) +
                   labels = c("-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5", "+6")) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.3)) +
  guides(fill = FALSE)
  

  
pos_colour <- rep("grey30", 7)
pos_colour[5:6] <- "red"
  

p_acceptor_count <- ggplot(splicing_mut_info_filt_snv_count %>% 
                             filter(Type_Motif == "acceptor", GenomonSplicingMutation != "Normal splicing"), 
                           aes(x = Rel_Start_Motif, y = count, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", y = "SAV count", fill = "") +
  ggtitle("Acceptor") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(# axis.text.x = element_text(colour = pos_colour),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:7, 
                   # labels = c("Y", "Y", "N", "C", "A", "G", "G")) +
                   labels = c("+6", "+5", "+4", "+3", "+2", "+1", "-1")) + 
  scale_y_continuous(expand = c(0, 0), limit = c(0, 3000)) +
  guides(fill = FALSE)


p_acceptor_ratio <- ggplot(splicing_mut_info_filt_snv_ratio %>% 
         filter(Type_Motif == "acceptor", GenomonSplicingMutation != "Normal splicing"), 
       aes(x = Rel_Start_Motif, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", y = "SAV ratio", fill = "") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(# axis.text.x = element_text(colour = pos_colour),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:7, 
                   # labels = c("Y", "Y", "N", "C", "A", "G", "G")) +
                   labels = c("+6", "+5", "+4", "+3", "+2", "+1", "-1")) + 
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.3)) + 
  guides(fill = FALSE)




g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


p_dummy_for_legend <- ggplot(splicing_mut_info_filt_snv_ratio %>%
                               filter(GenomonSplicingMutation != "Normal splicing"), 
                             aes(x = Rel_Start_Motif, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", fill = "") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() + 
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow=2, byrow=TRUE))


p_donor <- plot_grid(p_donor_count, p_donor_ratio, ncol = 1, rel_heights = c(1, 0.9), align = "v")
p_acceptor <- plot_grid(p_acceptor_count, p_acceptor_ratio, ncol = 1, rel_heights = c(1, 0.9), align = "v")

plot_grid(plot_grid(p_donor, p_acceptor, ncol = 2, align = "h", rel_widths = c(1, 0.9)), 
          ggdraw() + draw_label(" ", size = 7),
          g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(2, 0.15, 0.15))

# p_donor_acceptor_count <- plot_grid(p_donor_count, p_acceptor_count, ncol = 2, rel_widths = c(1, 0.9), align = "h")
# p_donor_acceptor_ratio <- plot_grid(p_donor_ratio, p_acceptor_ratio, ncol = 2, rel_widths = c(1, 0.9), align = "h")
# plot_grid(p_donor_acceptor_count, p_donor_acceptor_ratio, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.95, 0.2), align = "v")

write.table(splicing_mut_info_filt_snv_ratio, "../table/splicing_snv_ratio.txt", quote = FALSE, row.names = FALSE, sep = "\t")

ggsave("../figure/splicing_snv_ratio.tiff", width = 10, height = 7, dpi = 600, units = "cm")


##########
##########
# indel canonical or non-canonical

splicing_mut_info_filt_indel <- splicing_mut_info_filt %>% filter(Ref_Mut == "-" | Alt_Mut == "-")

Is_Canonical <- rep("non-canonical", nrow(splicing_mut_info_filt_indel))

# insertion
Is_Canonical[splicing_mut_info_filt_indel$Ref_Mut == "-" & 
               splicing_mut_info_filt_indel$Type_Motif == "donor" & 
               splicing_mut_info_filt_indel$Strand_Motif == "+" &
               splicing_mut_info_filt_indel$Rel_Start_Motif == 4] <- "canonical"

# Is_Canonical[splicing_mut_info_filt_indel$Ref_Mut == "-" & 
#                splicing_mut_info_filt_indel$Type_Motif == "donor" & 
#                splicing_mut_info_filt_indel$Strand_Motif == "-" &
#                splicing_mut_info_filt_indel$Rel_Start_Motif == 5] <- "canonical"

Is_Canonical[splicing_mut_info_filt_indel$Ref_Mut == "-" & 
               splicing_mut_info_filt_indel$Type_Motif == "acceptor" & 
               splicing_mut_info_filt_indel$Strand_Motif == "+" &
               splicing_mut_info_filt_indel$Rel_Start_Motif == 5] <- "canonical"

# Is_Canonical[splicing_mut_info_filt_indel$Ref_Mut == "-" & 
#                splicing_mut_info_filt_indel$Type_Motif == "acceptor" & 
#                splicing_mut_info_filt_indel$Strand_Motif == "-" &
#                splicing_mut_info_filt_indel$Rel_Start_Motif == 6] <- "canonical"

# deletion
Is_Canonical[splicing_mut_info_filt_indel$Alt_Mut == "-" & 
               splicing_mut_info_filt_indel$Type_Motif == "donor" & 
               splicing_mut_info_filt_indel$Rel_Start_Motif <= 5 &
               splicing_mut_info_filt_indel$Rel_End_Motif >= 4] <- "canonical"

Is_Canonical[splicing_mut_info_filt_indel$Alt_Mut == "-" & 
               splicing_mut_info_filt_indel$Type_Motif == "acceptor" & 
               splicing_mut_info_filt_indel$Rel_Start_Motif <= 6 &
               splicing_mut_info_filt_indel$Rel_End_Motif >= 5] <- "canonical"

splicing_mut_info_filt_indel$Is_Canonical <- Is_Canonical

Indel_Type <- rep("", nrow(splicing_mut_info_filt_indel))
# Indel_Type[splicing_mut_info_filt_indel$Alt_Mut == "-" & Is_Canonical == "canonical"] <- "Del (C)"
# Indel_Type[splicing_mut_info_filt_indel$Alt_Mut == "-" & Is_Canonical == "non-canonical"] <- "Del (N)"
# Indel_Type[splicing_mut_info_filt_indel$Ref_Mut == "-" & Is_Canonical == "canonical"] <- "Ins (C)"
# Indel_Type[splicing_mut_info_filt_indel$Ref_Mut == "-" & Is_Canonical == "non-canonical"] <- "Ins (N)"
Indel_Type[splicing_mut_info_filt_indel$Alt_Mut == "-" & Is_Canonical == "canonical"] <- "Canonical\ndeletion"
Indel_Type[splicing_mut_info_filt_indel$Alt_Mut == "-" & Is_Canonical == "non-canonical"] <- "Noncanonical\ndeletion"
Indel_Type[splicing_mut_info_filt_indel$Ref_Mut == "-" & Is_Canonical == "canonical"] <- "Canonical\ninsertion"
Indel_Type[splicing_mut_info_filt_indel$Ref_Mut == "-" & Is_Canonical == "non-canonical"] <- "Noncanonical\ninsertion"
 
Indel_Type <- factor(Indel_Type, levels = c("Canonical\ndeletion", "Noncanonical\ndeletion", "Canonical\ninsertion", "Noncanonical\ninsertion"))

splicing_mut_info_filt_indel$Indel_Type <- Indel_Type
  
# splicing_mut_info_filt_indel$InsDel <- rep("deletion", nrow(splicing_mut_info_filt_indel))
# splicing_mut_info_filt_indel$InsDel[splicing_mut_info_filt_indel$Ref_Mut == "-"] <- "insertion"


# obtain the ratio
splicing_mut_info_filt_indel_count <- splicing_mut_info_filt_indel %>% 
  group_by(Indel_Type, Type_Motif, GenomonSplicingMutation) %>% 
  summarize(count = n())

splicing_mut_info_filt_indel_count_total <- splicing_mut_info_filt_indel_count %>% 
  group_by(Indel_Type, Type_Motif) %>% 
  summarize(total_count = sum(count))

splicing_mut_info_filt_indel_ratio <- 
  left_join(splicing_mut_info_filt_indel_count, splicing_mut_info_filt_indel_count_total, by = c("Indel_Type", "Type_Motif")) %>%
  mutate(ratio = count / total_count) 


write.table(splicing_mut_info_filt_indel_ratio,
  "../temporary/splicing_indel_ratio.txt", quote = FALSE, row.names = FALSE, sep = "\t")


p_donor_indel_count <- ggplot(splicing_mut_info_filt_indel_count %>% 
         filter(Type_Motif == "donor", GenomonSplicingMutation != "Normal splicing"), 
       aes(x = Indel_Type, y = count, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  ggtitle("Donor") +
  labs(x = "", y = "SAV count") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(# axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # strip.text = element_text(size = rel(1.2))) +
        # legend.position = "bottom") +
        ) + 
  scale_y_continuous(expand = c(0, 0), limit = c(0, 300)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  guides(fill = FALSE)

p_donor_indel_ratio <- ggplot(splicing_mut_info_filt_indel_ratio %>% 
                                filter(Type_Motif == "donor", GenomonSplicingMutation != "Normal splicing"), 
                              aes(x = Indel_Type, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = splicing_class_colour) +
  labs(x = "", y = "SAV ratio") +
  my_theme() +
  theme(# axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # strip.text = element_text(size = rel(1.2))) +
        # legend.position = "bottom") +
        ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.35)) +
  guides(fill = FALSE)

p_acceptor_indel_count <- ggplot(splicing_mut_info_filt_indel_count %>% 
                                filter(Type_Motif == "acceptor", GenomonSplicingMutation != "Normal splicing"), 
                              aes(x = Indel_Type, y = count, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = splicing_class_colour) +
  ggtitle("Acceptor") +
  labs(x = "", y = "SAV count") +
  my_theme() +
  theme(# axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # strip.text = element_text(size = rel(1.2))) +
        # legend.position = "bottom") +
        ) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 300)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  guides(fill = FALSE)

p_acceptor_indel_ratio <- ggplot(splicing_mut_info_filt_indel_ratio %>% 
                                filter(Type_Motif == "acceptor", GenomonSplicingMutation != "Normal splicing"), 
                              aes(x = Indel_Type, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  labs(x = "", y = "SAV ratio") +
  theme(# axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        # strip.text = element_text(size = rel(1.2))) +
        # legend.position = "bottom") +
        ) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.35)) +
  # guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  guides(fill = FALSE)

p_dummy_for_legend2 <- ggplot(splicing_mut_info_filt_snv_ratio %>%
                               filter(GenomonSplicingMutation != "Normal splicing"),
                             aes(x = Rel_Start_Motif, y = ratio, fill = GenomonSplicingMutation)) +
  geom_bar(stat = "identity") +
  labs(x = "", fill = "") +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() 

# p_donor_acceptor_indel_count <- plot_grid(p_donor_indel_count, p_acceptor_indel_count, ncol = 2, align = "h")
# p_donor_acceptor_indel_ratio <- plot_grid(p_donor_indel_ratio, p_acceptor_indel_ratio, ncol = 2, align = "h")
# plot_grid(p_donor_acceptor_indel_count, p_donor_acceptor_indel_ratio, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 1, 0.2))

p_donor_indel <- plot_grid(p_donor_indel_count, p_donor_indel_ratio, ncol = 1, rel_heights = c(1, 0.9), align = "v")
p_acceptor_indel <- plot_grid(p_acceptor_indel_count, p_acceptor_indel_ratio, ncol = 1, rel_heights = c(1, 0.9), align = "v")

# plot_grid(plot_grid(p_donor_indel, p_acceptor_indel, ncol = 2, align = "h"), g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(2, 0.2))
plot_grid(plot_grid(p_donor_indel, p_acceptor_indel, ncol = 2, align = "h"), g_legend(p_dummy_for_legend2), ncol = 2, rel_widths = c(1, 0.3))


ggsave("../figure/splicing_indel_ratio.tiff", width = 12, height = 9, dpi = 600, units = "cm")



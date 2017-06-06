library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("../../../conf/plot_config.R")

GC_info <- read.table("../temporary/TCGA.savnet.gc_intron.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

GC_info <- GC_info %>% filter(Splice_Class != "no-change" | FPKM >= 10)

print(nrow(GC_info))

GC_info_proc <- GC_info %>% 
  select( Mutation_Key, Type_Motif, Splice_Class, GC_intron_5prime, GC_exon, GC_intron_3prime) %>%
  mutate(GC_intron = 0.5 * (GC_intron_5prime + GC_intron_3prime)) %>%
  mutate(GC_exon_intron_diff = GC_exon - GC_intron) %>%
  distinct() %>%
  gather(Is_Intron, GC_ratio, GC_exon, GC_intron, GC_exon_intron_diff, GC_intron_5prime, GC_intron_3prime)
  # gather(Pos, GC_ratio, -Type_Motif, -Splice_Class)
  

GC_info_proc$Type_Motif <- 
  factor(GC_info_proc$Type_Motif,
         levels = c("donor", "acceptor"),
         labels = c("Donor", "Acceptor"))

  
# GC_info_proc$Splice_Class[
#   GC_info_proc$Splice_Class %in% c("intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site", 
#                                    "alternative-5'-splice-site", "alternative-3'-splice-site")] <- "alternative-splice-site"

GC_info_proc$Splice_Class[
  GC_info_proc$Splice_Class %in% c("alternative-5'-splice-site", "intronic-alternative-5'-splice-site")] <- 
    "alternative-5'-splice-site"

GC_info_proc$Splice_Class[
  GC_info_proc$Splice_Class %in% c("alternative-3'-splice-site", "intronic-alternative-3'-splice-site")] <- 
    "alternative-3'-splice-site"


GC_info_proc$Splice_Class[GC_info_proc$Splice_Class == "opposite-side-intron-retention"] <- "intron-retention"


GC_info_proc$Splice_Class2 <- factor(GC_info_proc$Splice_Class,
                                    levels = c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site", 
                                               "intron-retention", "complex", "no-change"),
                                    labels = c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                                               "Intron retention", "Complex", "Normal splicing"))

GC_info_proc$Is_Intron2 <- factor(GC_info_proc$Is_Intron, 
                                 levels = c("GC_intron_5prime", "GC_exon", "GC_intron_3prime", "GC_intron", "GC_exon_intron_diff"),
                                 labels = c("5' intron", "Exon", "3' intron", "Intron", "Exon intron diff"))



g_gc_d <- ggplot(GC_info_proc %>% 
  filter(Type_Motif == "Donor", Splice_Class2 != "Alternative 3'SS", Is_Intron2 %in% c("5' intron", "Exon", "3' intron")), 
  aes(x = Is_Intron2, y = GC_ratio, fill = Splice_Class2)) + 
  geom_boxplot(outlier.size = 0.3, size = 0.3) +
  labs(x = "", y = "GC content") +
  ggtitle("Donor disruption") +
  facet_grid(.~Splice_Class2) +
  guides(fill = FALSE) +
  ylim(c(0.2, 0.8)) + 
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 5)) 

g_gc_a <- ggplot(GC_info_proc %>% 
  filter(Type_Motif == "Acceptor", Splice_Class2 != "Alternative 5'SS", Is_Intron2 %in% c("5' intron", "Exon", "3' intron")),
  aes(x = Is_Intron2, y = GC_ratio, fill = Splice_Class2)) + 
  geom_boxplot(outlier.size = 0.3, size = 0.3) +
  labs(x = "", y = "GC content") +
  ggtitle("Acceptor disruption") +
  facet_grid(.~Splice_Class2) +
  guides(fill = FALSE) +
  ylim(c(0.2, 0.8)) + 
  scale_fill_manual(values = splicing_class_colour) + 
  my_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 5))

plot_grid(g_gc_d, g_gc_a, ncol = 1)

ggsave("../figure/diff_gc_content.tiff", width = 9, height = 11, dpi = 600, units = "cm")




# GC_info_proc$Splice_Class3 <- factor(GC_info_proc$Splice_Class,
#                                      levels = rev(c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
#                                                 "intron-retention", "complex", "no-change")),
#                                      labels = rev(c("Exon skipping", "Alternative 5' splice site", "Alternative 3' splice site",
#                                                 "Intron retention", "Complex", "Normal splicing")))


g_dgc_d <- ggplot(GC_info_proc %>% 
  filter(Is_Intron == "GC_exon_intron_diff" & Type_Motif == "Donor", Splice_Class2 != "Alternative 3'SS"),
  aes(x = Splice_Class2, y = GC_ratio, fill = Splice_Class2)) + 
  geom_boxplot(outlier.size = 0.3, size = 0.3) +
  labs(x = "", y = "Diff. of GC content between exonic and intronic regions") +
  ggtitle("Donor disruption") +
  # facet_grid(Type_Motif~.) +
  guides(fill = FALSE) +
  coord_flip() + 
  my_theme() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
  #       axis.title = element_text(size = 6),
  #       strip.text = element_text(size = 5)) +
  scale_fill_manual(values = splicing_class_colour) 

g_dgc_a <- ggplot(GC_info_proc %>% 
  filter(Is_Intron == "GC_exon_intron_diff" & Type_Motif == "Acceptor", Splice_Class2 != "Alternative 5'SS"),
  aes(x = Splice_Class2, y = GC_ratio, fill = Splice_Class2)) + 
  geom_boxplot(outlier.size = 0.3, size = 0.3) +
  labs(x = "", y = "Diff. of GC content between exonic and intronic regions") +
  ggtitle("Acceptor disruption") +
  # facet_grid(Type_Motif~.) +
  guides(fill = FALSE) +
  coord_flip() + 
  my_theme() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
  #       axis.title = element_text(size = 6),
  #       strip.text = element_text(size = 5)) +
  scale_fill_manual(values = splicing_class_colour)

plot_grid(g_dgc_d, g_dgc_a, ncol = 1)


ggsave("../figure/diff_exon_intron_diff_gc_content.tiff", width = 12, height = 7, dpi = 600, units = "cm")



##########

Len_info_proc <- GC_info %>%
  select(Mutation_Key, Type_Motif, Splice_Class, Len_intron_5prime, Len_exon, Len_intron_3prime) %>%
  distinct() %>%
  gather(key = Is_Intron, value = Len, Len_intron_5prime, Len_exon, Len_intron_3prime)

Len_info_proc$Type_Motif <-
  factor(Len_info_proc$Type_Motif,
         levels = c("donor", "acceptor"),
         labels = c("Donor", "Acceptor"))

Len_info_proc$Splice_Class[
  Len_info_proc$Splice_Class %in% c("alternative-5'-splice-site", "intronic-alternative-5'-splice-site")] <-
    "alternative-5'-splice-site"

Len_info_proc$Splice_Class[
  Len_info_proc$Splice_Class %in% c("alternative-3'-splice-site", "intronic-alternative-3'-splice-site")] <-
    "alternative-3'-splice-site"


Len_info_proc$Splice_Class[Len_info_proc$Splice_Class == "opposite-side-intron-retention"] <- "intron-retention"


Len_info_proc$Splice_Class2 <- factor(Len_info_proc$Splice_Class,
                                    levels = c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
                                               "intron-retention", "complex", "no-change"),
                                    labels = c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                                               "Intron retention", "Complex", "Normal splicing"))

Len_info_proc$Is_Intron2 <- factor(Len_info_proc$Is_Intron,
                                 levels = c("Len_intron_5prime", "Len_exon", "Len_intron_3prime"),
                                 labels = c("5' intron", "Exon", "3' intron"))


g_len_d <- ggplot(Len_info_proc %>%
  filter(Type_Motif == "Donor", Splice_Class2 != "Alternative 3'SS", Is_Intron2 %in% c("5' intron", "Exon", "3' intron")),
  aes(x = Is_Intron2, y = log10(Len), fill = Splice_Class2)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3) +
  labs(x = "", y = "Log10(length)") +
  ggtitle("Donor disruption") +
  facet_grid(.~Splice_Class2) +
  guides(fill = FALSE) +
  ylim(c(1.5, 5)) +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 5))

g_len_a <- ggplot(Len_info_proc %>%
  filter(Type_Motif == "Acceptor", Splice_Class2 != "Alternative 5'SS", Is_Intron2 %in% c("5' intron", "Exon", "3' intron")),
  aes(x = Is_Intron2, y = log10(Len), fill = Splice_Class2)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3) +
  labs(x = "", y = "Log10(length)") +
  ggtitle("Acceptor disruption") +
  facet_grid(.~Splice_Class2) +
  guides(fill = FALSE) +
  ylim(c(1.5, 5)) +
  scale_fill_manual(values = splicing_class_colour) +
  my_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 5))

plot_grid(g_len_d, g_len_a, ncol = 2)

ggsave("../figure/diff_length.tiff", width = 18, height = 6, dpi = 600, units = "cm")


library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

source("../../../conf/plot_config.R")

##########
# gene summary
splicing_mutation <- read.table("../../output/rescue/TCGA.savnet.with_rescued.result.txt", 
    header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE) %>% 
    filter(IR_filtered != "TRUE")

# splicing_mutation <- read.table("omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)


# cancer_gene
cancer_gene <- splicing_mutation %>% 
  filter(Is_Cancer_Gene == "TRUE") %>% 
  select(Gene_Symbol) %>% 
  distinct()

splicing_mutation_count <- splicing_mutation %>% 
  # filter(Is_Cancer_Gene == "TRUE") %>% 
  select(Cancer_Type, Sample_Name, Gene_Symbol) %>%
  distinct() %>%
  group_by(Cancer_Type, Gene_Symbol) %>%
  summarize(count = n()) %>% arrange(desc(count))

splicing_mutation_count_total <- splicing_mutation_count %>%
  group_by(Gene_Symbol) %>%
  summarize(total_count = sum(count)) %>%
  filter(total_count >= 10) %>% 
  arrange(total_count)


splicing_mutation_count_proc <- splicing_mutation_count %>% 
  filter(Gene_Symbol %in% splicing_mutation_count_total$Gene_Symbol)

splicing_mutation_count_proc$Gene_Symbol2 <-
  factor(splicing_mutation_count_proc$Gene_Symbol,
         levels = splicing_mutation_count_total$Gene_Symbol)

gene_colour <- rep("grey30", nrow(splicing_mutation_count_proc))
gene_colour[levels(splicing_mutation_count_proc$Gene_Symbol2) %in% cancer_gene$Gene_Symbol] <- "#b2182b"


p_summary <- ggplot(splicing_mutation_count_proc,
       aes(x = Cancer_Type, y = Gene_Symbol2, size = count)) + 
  geom_point(colour = "#2166ac") +
  labs(size = "Number of cases") +
  theme_minimal() +
  ggtitle("Frequency per cancer type") +
  labs(x = "Cancer Type", y = "Gene") +
  scale_size(range = c(1, 4)) +
  theme(title = element_text(size = 7),
        axis.text = element_text(size = 7), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(colour =  gene_colour, face = "italic"),
        legend.position = "bottom")


p_total_count <- ggplot(splicing_mutation_count_total, aes(x = reorder(Gene_Symbol, total_count), y = total_count)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + theme_minimal() +
  ggtitle("Total frequency") +
  labs(x = "", y = "") +
  theme(axis.text.y = element_text(colour =  gene_colour),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

##########
# sasm type info

sasm_type_info <- splicing_mutation

sasm_type_info[grep("creation", sasm_type_info$Mutation_Type), 
  "Mutation_Type"] <- "motif_creation"

sasm_type_info[
  intersect(grep("disruption", sasm_type_info$Mutation_Type),
            which(sasm_type_info$Is_Canonical == "canonical")),
  "Mutation_Type"] <- "motif_disruption_C"

sasm_type_info[
  intersect(grep("disruption", sasm_type_info$Mutation_Type),
            which(sasm_type_info$Is_Canonical == "non-canonical")),
  "Mutation_Type"] <- "motif_disruption_N"



sasm_type_info_proc <- sasm_type_info %>% 
  group_by(Sample_Name, Gene_Symbol, Mutation_Key, Mutation_Type) %>% 
  distinct() %>% 
  group_by(Gene_Symbol, Mutation_Type) %>% 
  summarize(count = n()) %>% 
  spread(key = Mutation_Type, value = count, fill = 0) %>%
  mutate(total_count = motif_creation + motif_disruption_C + motif_disruption_N) %>%
  mutate(motif_creation = motif_creation / total_count,
         motif_disruption_C = motif_disruption_C / total_count,
         motif_disruption_N = motif_disruption_N / total_count) %>%
  gather(key = mut_type, value = ratio, motif_creation, motif_disruption_C, motif_disruption_N) %>%
  filter(Gene_Symbol %in% splicing_mutation_count_total$Gene_Symbol)

sasm_type_info_proc$Gene_Symbol2 <- factor(sasm_type_info_proc$Gene_Symbol,
                                          levels = splicing_mutation_count_total$Gene_Symbol)


sasm_type_info_proc$mut_type <- factor(sasm_type_info_proc$mut_type,
                                      levels = rev(c("motif_disruption_C", "motif_disruption_N", "motif_creation")),
                                      labels = rev(c("Canonical site disruption", "Noncanonical site disruption", "Creation")))


p_sasmtype <- ggplot(sasm_type_info_proc, aes(x = Gene_Symbol2, y = ratio, fill = mut_type)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Mutation type") +
  coord_flip() + theme_minimal() +
  labs(x = "", y = "Relative frequency", fill = "Mutation type") +
  theme_minimal() +
  theme(
    title = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6),
    legend.key.size = unit(0.25, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(colour =  gene_colour, face = "italic"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom") +
  scale_fill_manual(breaks = c("Canonical site disruption", "Noncanonical site disruption", "Creation"),
                    values = c("Canonical site disruption" = "#a6d854", "Noncanonical site disruption" = "#ffd92f", "Creation" = "#e78ac3"))



##########
# splicing class

sp_class_info <- splicing_mutation %>% 
  group_by(Sample_Name, Gene_Symbol, Mutation_Key) %>% 
  summarize(Splicing_Class_Paste = paste(Splicing_Class, collapse = ";"),
            Supporting_Read_Num_Paste = paste(Supporting_Read_Num, collapse = ","))
            

splcing_class2 <- sp_class_info$Splicing_Class_Paste

splcing_class2 <- unlist(
  lapply(
    1:nrow(sp_class_info),
    function(x) {
      treads <- as.numeric(strsplit(sp_class_info$Supporting_Read_Num_Paste[x], ",")[[1]])
      tind <- which(treads == max(treads))[1]
      return(strsplit(sp_class_info$Splicing_Class_Paste[x], ";")[[1]][tind])
    }
  )
)

splcing_class2[splcing_class2 == "opposite-side-intron-retention"] <- "intron-retention"
splcing_class2[splcing_class2 == "intronic-alternative-5'-splice-site"] <- "alternative-5'-splice-site"
splcing_class2[splcing_class2 == "intronic-alternative-3'-splice-site"] <- "alternative-3'-splice-site"

sp_class_info$Splicing_Class2 <- splcing_class2

sp_class_info_proc <- sp_class_info %>%
  group_by(Gene_Symbol, Splicing_Class2) %>% summarize(count = n()) %>%
  spread(key = Splicing_Class2, value = count, fill = 0) %>%
  mutate(total_count = `exon-skip` + `alternative-5'-splice-site` + `alternative-3'-splice-site` +
           `intron-retention`) %>%
  mutate(`exon-skip` = `exon-skip` / total_count,
         `alternative-5'-splice-site` = `alternative-5'-splice-site` / total_count,
         `alternative-3'-splice-site` = `alternative-3'-splice-site` / total_count,
         `intron-retention` = `intron-retention` / total_count) %>%
  gather(key = splice_class, value = ratio,
         `exon-skip`, `alternative-5'-splice-site`, `alternative-3'-splice-site`, 
         `intron-retention`) %>%
  filter(Gene_Symbol %in% splicing_mutation_count_total$Gene_Symbol)

sp_class_info_proc$Gene_Symbol2 <- factor(sp_class_info_proc$Gene_Symbol,
                                         levels = splicing_mutation_count_total$Gene_Symbol)

sp_class_info_proc$splice_class2 <- factor(sp_class_info_proc$splice_class,
                                          levels = rev(c("exon-skip", "alternative-5'-splice-site",
                                                     "alternative-3'-splice-site", "intron-retention")),
                                          labels = rev(c("Exon skip", "Alternative 5'-ss",
                                                     "Alternative 3'-ss", "Intron retention")))


p_spliceclass <- ggplot(sp_class_info_proc, aes(x = Gene_Symbol2, y = ratio, fill = splice_class2)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + theme_minimal() +
  ggtitle("Splcing class") +
  labs(x = "", y = "Relative frequency", fill = "Splicing class") +
  theme_minimal() +
  theme(
    title = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6),
    legend.key.size = unit(0.25, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(colour =  gene_colour, face = "italic"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom") +
  scale_fill_manual(values = splicing_class_colour, 
                    breaks = c("Exon skip", "Alternative 5'-ss",
                              "Alternative 3'-ss", "Intron retention"))


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p_summary_legend <- g_legend(p_summary)
p_sasmtype_legend <- g_legend(p_sasmtype)
p_spliceclass_legend <- g_legend(p_spliceclass)

p_main_panel <- plot_grid(p_summary + guides(size = FALSE), 
          p_sasmtype + guides(fill = FALSE), 
          p_spliceclass + guides(fill = FALSE),
          nrow = 1, rel_widths = c(0.6, 0.2, 0.2),
          align = "h")

# p_legend_panel <- plot_grid(
#   plot_grid(p_summary_legend, p_sasmtype_legend, nrow = 1, align = "h"),
#   p_spliceclass_legend, nrow = 2, align = "h")

p_legend_panel <- plot_grid(p_summary_legend, p_sasmtype_legend, p_spliceclass_legend, nrow = 3)
  
plot_grid(p_main_panel, p_legend_panel, nrow = 2,
          rel_heights = c(0.90, 0.1))

  
ggsave("../figure/cancer_gene_summary.tiff", width = 20, height = 20, dpi = 600, units = "cm")





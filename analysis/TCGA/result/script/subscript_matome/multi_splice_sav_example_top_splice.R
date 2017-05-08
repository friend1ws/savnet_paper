library(ggplot2)
library(dplyr)
library(cowplot)

source("subscript_matome/sav_profile_function.R")

refGene <- read.table("../../../db/refGene/refGene.txt.gz",
                      sep = "\t", header = FALSE, stringsAsFactors = FALSE)

##########

splicing_mutation <- read.table("../temporary/TCGA.motif_read_num.txt", header = TRUE, sep = "\t", 
                                quote="", stringsAsFactors = FALSE) %>% filter(Supporting_Read_Num != 0)

splicing_mutation[splicing_mutation$Splicing_Class == "intronic-alternative-5'-splice-site", "Splicing_Class"] <- "alternative-5'-splice-site"
splicing_mutation[splicing_mutation$Splicing_Class == "intronic-alternative-3'-splice-site", "Splicing_Class"] <- "alternative-3'-splice-site"
splicing_mutation[splicing_mutation$Splicing_Class == "opposite-side-intron-retention", "Splicing_Class"] <- "intron-retention"

splicing_mutation$Splicing_Class <-
  factor(splicing_mutation$Splicing_Class,
         levels = c("exon-skip", "alternative-5'-splice-site",
                    "alternative-3'-splice-site",
                    "intron-retention"),
         labels = c("Exon skipping", "Alternative 5'SS",
                    "Alternative 3'SS", "Intron retention"))

splicing_mutation$Is_Inframe <-
  factor(splicing_mutation$Is_Inframe,
         levels = c("in-frame", "---"),
         labels = c("In-frame", "Frameshift"))


p_TP53 <- get_print_info(splicing_mutation, "TCGA-85-8479-01", "TP53", "17,7579312,C,A", "17:7578554-7579311") + guides(color = FALSE, linetype = FALSE)
p_GATA3 <- get_print_info(splicing_mutation, "TCGA-B6-A0WS-01", "GATA3", "10,8111432,TCA,T", "10:8106102-8111442") + guides(color = FALSE, linetype = FALSE)



label1 <- ggdraw() + draw_label(" ", size = 6)
label2 <- ggdraw() + draw_label("Fraction of the most frequent\ntype of abnormal splicing", size = 6)
label3 <- ggdraw() + 
            draw_label(get_most_freq_splicing_calc_expression(splicing_mutation, "TCGA-85-8479-01", "TP53", "17,7579312,C,A", "17:7578554-7579311"), size = 6)
label4 <- ggdraw() + 
            draw_label(get_most_freq_splicing_calc_expression(splicing_mutation, "TCGA-B6-A0WS-01", "GATA3", "10,8111432,TCA,T", "10:8106102-8111442"), size = 6)

p_legend1 <- g_legend(get_legend_info() + theme(legend.justification = 0.1))
p_legend2 <- g_legend(ggplot(data.frame(x = c(1), y = c(1), col = c("The most frequent type of abnormal splicing")), 
                        aes(x = x, y = y, color = col)) +
                      geom_point(shape = 19) +
                      theme(legend.position = "bottom",
                            legend.text = element_text(size = 6),
                            legend.key = element_blank(),
                            legend.key.size = unit(0.25, "cm"),
                            legend.justification = 0.045,
                            legend.margin = margin(0.5, 0.5, 0.5, 0.5)) +
                      scale_colour_manual(values = c("The most frequent type of abnormal splicing" = "#ff7f00")) +
                      labs(color = ""))
 


plot_grid(plot_grid(label1, p_TP53, p_GATA3, align = "h", nrow = 1, rel_widths = c(0.4, 1, 1)),
          plot_grid(label2, label3, label4, nrow = 1, rel_widths = c(0.4, 1, 1)),
          plot_grid(label1, p_legend1, nrow = 1, rel_widths = c(0.4, 2)),
          plot_grid(label1, p_legend2, nrow = 1, rel_widths = c(0.4, 2)), ncol = 1, rel_heights = c(0.8, 0.2, 0.12, 0.12))


ggsave("../figure/multi_splice_sav_example_top_splice.tiff", width = 20, height = 6, dpi = 600, units = "cm")


plot_grid(p_legend1)
ggsave("../figure/profile_legend.tiff", width = 16, height = 0.8, dpi = 600, units = "cm")


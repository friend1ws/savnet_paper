library(ggplot2)
library(dplyr)
library(cowplot)

# source("../../../conf/plot_config.R")
source("subscript_matome/sav_profile_function.R")

refGene <- read.table("../../../db/refGene/refGene.txt.gz",
                      sep = "\t", header = FALSE, stringsAsFactors = FALSE)

splicing_mutation <- read.table("../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt", header = TRUE, sep = "\t", 
                                quote="", stringsAsFactors = FALSE)
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


p_EIF1 <- get_print_info(splicing_mutation, "TCGA-A3-3313-01", "EIF1", "17,39846195,T,C") + guides(color = FALSE, linetype = FALSE)
p_MEF2B <- get_print_info(splicing_mutation, "TCGA-G8-6914-01", "MEF2B", "19,19257443,ACTGTAGAGGCTTCTCTGTG,A") + guides(color = FALSE, linetype = FALSE)
p_POLD2 <- get_print_info(splicing_mutation, "TCGA-AD-A5EJ-01", "POLD2", "7,44155494,T,C") + guides(color = FALSE, linetype = FALSE)
p_CCNG2 <- get_print_info(splicing_mutation, "TCGA-D8-A140-01", "CCNG2", "4,78082129,G,T") + guides(color = FALSE, linetype = FALSE)

# p_legend <- g_legend(get_print_info(splicing_mutation, "TCGA-AD-A5EJ-01", "POLD2", "7,44155494,T,C") +
#                guides(color = guide_legend(order = 1), linetype = guide_legend(order = 0)))
p_legend <- g_legend(get_legend_info())


plot_grid(plot_grid(p_MEF2B, p_EIF1, p_POLD2, p_CCNG2, ncol = 2),
          p_legend,
          # g_legend(get_print_info("POLD2", "7,44155494,T,C")),
          ncol = 1, rel_heights = c(0.9, 0.1))

ggsave("../figure/multi_splice_mutation_example.tiff", width = 20, height = 8, dpi = 600, units = "cm")


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


print("OK")

label1 <- ggdraw() + draw_label(" ", size = 6)
label2 <- ggdraw() + draw_label("The fraction of most frequent\ntype of abnormal splicing", size = 6)
label3 <- ggdraw() + 
            draw_label(get_most_freq_splicing_calc_expression(splicing_mutation, "TCGA-85-8479-01", "TP53", "17,7579312,C,A", "17:7578554-7579311"), size = 6)
label4 <- ggdraw() + 
            draw_label(get_most_freq_splicing_calc_expression(splicing_mutation, "TCGA-B6-A0WS-01", "GATA3", "10,8111432,TCA,T", "10:8106102-8111442"), size = 6)

p_legend <- g_legend(get_legend_info())

plot_grid(plot_grid(label1, p_TP53, p_GATA3, align = "h", nrow = 1, rel_widths = c(0.4, 1, 1)),
          plot_grid(label2, label3, label4, nrow = 1, rel_widths = c(0.4, 1, 1)),
          p_legend, ncol = 1, rel_heights = c(0.8, 0.2, 0.2))


ggsave("../figure/multi_splice_sav_example_top_splice.tiff", width = 20, height = 6, dpi = 600, units = "cm")


plot_grid(p_legend)
ggsave("../figure/profile_legend.tiff", width = 16, height = 0.8, dpi = 600, units = "cm")


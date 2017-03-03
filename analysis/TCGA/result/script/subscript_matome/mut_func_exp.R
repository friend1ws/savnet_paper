library(dplyr)
library(ggplot2)

source("../../../conf/plot_config.R")

mut_func_exp <- read.table("../../output/mut_func_exp/TCGA.mut_func_exp.merged.txt", sep = "\t", header = TRUE, 
                as.is=TRUE, quote="", stringsAsFactors = FALSE)

cancer_gene <- read.table("../../../db/cancer_gene/cancer_gene.txt", sep = "\t", header = TRUE) %>% 
    filter(Cancer_Gene_Census != "---")
cancer_gene_list <- as.character(cancer_gene[,"Gene_Symbol"])

print(cancer_gene_list)

func_class2 <- rep("other", nrow(mut_func_exp))

func_class2[
  mut_func_exp$Func_Class %in% c("frameshift deletion", "frameshift insertion")] <- 
  "Frameshift indel"

func_class2[
  mut_func_exp$Func_Class %in% c("nonframeshift deletion", "nonframeshift insertion")] <- 
  "Inframe indel"

func_class2[
  mut_func_exp$Func_Class %in% c("exon-skip")] <-
  "Exon skip (frameshift)"

func_class2[
  func_class2 == "Exon skip (frameshift)" &
    mut_func_exp$Is_Inframe == "in-frame"] <-
  "Exon skip (inframe)"

func_class2[
  mut_func_exp$Func_Class %in% c("alternative-5'-splice-site", "intronic-alternative-5'-splice-site")] <-
  "Alternative 5'-ss (frameshift)"

func_class2[
  mut_func_exp$Func_Class %in% c("alternative-3'-splice-site", "intronic-alternative-3'-splice-site")] <-
  "Alternative 3'-ss (frameshift)"

func_class2[
  func_class2 == "Alternative 5'-ss (frameshift)" &
  mut_func_exp$Is_Inframe == "in-frame"] <-
  "Alternative 5'-ss (inframe)"

func_class2[
  func_class2 == "Alternative 3'-ss (frameshift)" &
  mut_func_exp$Is_Inframe == "in-frame"] <-
  "Alternative 3'-ss (inframe)"


func_class2[
  mut_func_exp$Func_Class %in% c("intron-retention", "opposite-side-intron-retention")] <-
  "Intron retention"

func_class2[
  mut_func_exp$Func_Class %in% c("synonymous SNV")] <-
  "Silent"

func_class2[
  mut_func_exp$Func_Class %in% c("nonsynonymous SNV")] <-
  "Missense"

func_class2[
  mut_func_exp$Func_Class %in% c("stopgain")] <-
  "Nonsense"

func_class2[
  mut_func_exp$Func_Class %in% c("complex")] <- "Complex"

mut_func_exp$Func_Class2 <- factor(func_class2, 
                                   levels =rev(c("Silent", "Missense", "Nonsense", "Inframe indel", "Frameshift indel",
                                                 "Exon skip (frameshift)", "Exon skip (inframe)",
                                                 "Alternative 5'-ss (frameshift)", "Alternative 5'-ss (inframe)",
                                                 "Alternative 3'-ss (frameshift)", "Alternative 3'-ss (inframe)",
                                                 "Intron retention", "Complex", "other")))


 

ggplot(mut_func_exp %>% 
         filter(!Func_Class2 %in% c("other", "Inframe indel", "Frameshift indel"),
                FPKM_mean >= 10.0),
                # Gene_Symbol %in% mut_func_exp_gene_count$Gene_Symbol),
       aes(x = Func_Class2, y = FPKM_normalized, fill = Func_Class2)) + 
  geom_boxplot(outlier.size = 0.6, size = 0.3) + ylim(c(-3, 3)) + 
  coord_flip() + 
  labs(x = "", y = "Normalized FPKM") +
  my_theme() +
  scale_fill_manual(values = splicing_class_colour) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = FALSE)
  

ggsave("../figure/mut_fun_exp_fpkm10.pdf", width = 6, height = 2.8)



ggplot(mut_func_exp %>% 
         filter(!Func_Class2 %in% c("other", "Inframe indel", "Frameshift indel"),
                FPKM_mean >= 0.0),
       # Gene_Symbol %in% mut_func_exp_gene_count$Gene_Symbol),
       aes(x = Func_Class2, y = FPKM_normalized, fill = Func_Class2)) + 
  geom_boxplot(outlier.size = 0.6, size = 0.3) + ylim(c(-3, 3)) + 
  coord_flip() + 
  labs(x = "", y = "Normalized FPKM") +
  my_theme() +
  scale_fill_manual(values = splicing_class_colour) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = FALSE)

ggsave("../figure/mut_fun_exp.pdf", width = 6, height = 2.8)



ggplot(mut_func_exp %>% 
         filter(!Func_Class2 %in% c("other", "Inframe indel", "Frameshift indel"),
                Gene_Symbol %in% cancer_gene_list),
       # Gene_Symbol %in% mut_func_exp_gene_count$Gene_Symbol),
       aes(x = Func_Class2, y = FPKM_normalized, fill = Func_Class2)) + 
  geom_boxplot(outlier.size = 0.6, size = 0.3) + ylim(c(-3, 3)) + 
  coord_flip() + 
  labs(x = "", y = "Normalized FPKM") +
  my_theme() +
  scale_fill_manual(values = splicing_class_colour) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = FALSE)



ggsave("../figure/mut_fun_exp_CG.pdf", width = 6, height = 2.8)




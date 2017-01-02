library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

cancer_gene <- read.table("/home/yshira/mysoftware/sv_utils/cancer_gene/cancer_gene.txt", sep = "\t", stringsAsFactors = FALSE) %>% 
  filter(V2 != "---" | V3 != "---") %>%
  filter(V2 != "---" | V3 != "T") %>% 
  select(gene = V1)

gene_list <- read.table("/home/yshira/mysoftware/sv_utils/resource/refGene.txt.gz", sep="\t", stringsAsFactors = FALSE) %>% 
  filter(substring(V2, 0, 3) == "NM_") %>%
  select(gene = V13) %>% distinct(gene) %>% mutate(cancer_gene = ifelse(gene %in% cancer_gene$gene, "CG", "nonCG")) 

mut_info <- read.table("../matome/omega.mutation.merged.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  filter(!str_detect(Gene.refGene, "-") & !str_detect(Gene.refGene, ",")) %>%
  full_join(gene_list, by = c("Gene.refGene" = "gene")) %>%
  filter(!is.na(cancer_gene))


syn_info <- mut_info %>% 
  filter(Func.refGene == "exonic" & ExonicFunc.refGene == "synonymous SNV") %>% 
  # select(Sample_Name, Gene.refGene, Func.refGene, ExonicFunc.refGene, cancer_gene) %>% distinct() %>%
  group_by(cancer_gene) %>% summarize(count = n()) %>%
  spread(key = cancer_gene, value = count) %>% mutate(ratio = CG / (CG + nonCG))
  

nonsyn_info <- mut_info %>% 
  filter(Func.refGene == "exonic" & ExonicFunc.refGene == "nonsynonymous SNV") %>% 
  # select(Sample_Name, Gene.refGene, Func.refGene, ExonicFunc.refGene, cancer_gene) %>% distinct() %>%
  group_by(cancer_gene) %>% summarize(count = n()) %>%
  spread(key = cancer_gene, value = count) %>% mutate(ratio = CG / (CG + nonCG))

none_info <- mut_info %>% 
  filter(Func.refGene == "exonic" & ExonicFunc.refGene == "stopgain") %>% 
  # select(Sample_Name, Gene.refGene, Func.refGene, ExonicFunc.refGene, cancer_gene) %>% distinct() %>%
  group_by(cancer_gene) %>% summarize(count = n()) %>%
  spread(key = cancer_gene, value = count) %>% mutate(ratio = CG / (CG + nonCG))

lof_info <- mut_info %>% 
  filter(Func.refGene == "exonic" & ExonicFunc.refGene %in% c("stopgain", "frameshift deletion", "frameshift insertion")) %>% 
  # select(Sample_Name, Gene.refGene, Func.refGene, ExonicFunc.refGene, cancer_gene) %>% distinct() %>%
  group_by(cancer_gene) %>% summarize(count = n()) %>%
  spread(key = cancer_gene, value = count) %>% mutate(ratio = CG / (CG + nonCG))

sp_info <- mut_info %>% 
  filter(str_detect(Func.refGene, "splicing")) %>% 
  # select(Sample_Name, Gene.refGene, Func.refGene, ExonicFunc.refGene, cancer_gene) %>% distinct() %>%
  group_by(cancer_gene) %>% summarize(count = n()) %>%
  spread(key = cancer_gene, value = count) %>% mutate(ratio = CG / (CG + nonCG))


gsm_info <- mut_info %>% 
  filter(GSM != "no-change") %>% 
  group_by(GSM, cancer_gene) %>% summarize(count = n()) %>%
  spread(key = cancer_gene, value = count) %>% mutate(ratio = CG / (CG + nonCG))


type_order <- c("silent", "missense", "loss-of-function",
                "exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site", 
                "intron-retention", "complex")

cg_ratio <- data.frame(
  type = factor(type_order, levels = rev(type_order)),
  ratio = c(syn_info$ratio, nonsyn_info$ratio, lof_info$ratio,
            unlist(lapply(type_order[4:8], function(x) {as.numeric(gsm_info[gsm_info$GSM == x, "ratio"])}))),
  is_gsm = c(rep("non_gsm", 3), rep("gsm", 5)))


ggplot(cg_ratio, aes(x = type, y = ratio, fill = is_gsm)) + 
  geom_bar(stat = "identity")  + 
  coord_flip() +
  labs(x = "", y = "Cancer gene ratio") +
  theme_minimal() +
  scale_fill_manual(values = c(gsm = "#33a02c", non_gsm = "#b2df8a")) +
  guides(fill = FALSE)


ggsave("../matome/cancer_gene_ratio.png", width = 6, height = 4)


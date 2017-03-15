library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

source("../../../conf/plot_config.R")

cancer_gene <- read.table("../../../db/cancer_gene/cancer_gene.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
  

gene_list <- read.table("../../../db/refGene/refGene.txt.gz", sep="\t", stringsAsFactors = FALSE) %>% 
  filter(substring(V2, 0, 3) == "NM_") %>%
  select(Gene_Symbol = V13) %>% distinct(Gene_Symbol) %>% left_join(cancer_gene, key = "Gene_Symbol") 


gene_list$LawrenceEtAl_2014[is.na(gene_list$LawrenceEtAl_2014)] <- "NonCG"
gene_list$LawrenceEtAl_2014[gene_list$LawrenceEtAl_2014 == "---"] <- "NonCG"
gene_list$LawrenceEtAl_2014[gene_list$LawrenceEtAl_2014 != "NonCG"] <- "CG"

gene_list$Cancer_Gene_Census[is.na(gene_list$Cancer_Gene_Census)] <- "NonCG"
gene_list$Cancer_Gene_Census[gene_list$Cancer_Gene_Census == "---"] <- "NonCG"
gene_list$Cancer_Gene_Census[gene_list$Cancer_Gene_Census != "NonCG"] <- "CG"

gene_list$VogelsteinEtAl_2013[is.na(gene_list$VogelsteinEtAl_2013)] <- "NonCG"
gene_list$VogelsteinEtAl_2013[gene_list$VogelsteinEtAl_2013 == "---"] <- "NonCG"
gene_list$VogelsteinEtAl_2013[gene_list$VogelsteinEtAl_2013 != "NonCG"] <- "CG"

gene_list$YeEtAl_2016[is.na(gene_list$YeEtAl_2016)] <- "NonCG"
gene_list$YeEtAl_2016[gene_list$YeEtAl_2016 == "---"] <- "NonCG"
gene_list$YeEtAl_2016[gene_list$YeEtAl_2016 != "NonCG"] <- "CG"



mut_info <- read.table("../temporary/TCGA.mutation.merged.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  filter(!str_detect(Gene.refGene, "-") & !str_detect(Gene.refGene, ",")) %>%
  inner_join(gene_list, by = c("Gene.refGene" = "Gene_Symbol")) 



mut_func <- rep("Other", nrow(mut_info))
mut_func[mut_info$Func.refGene == "exonic" & mut_info$ExonicFunc.refGene == "synonymous SNV"] <- "Silent"
mut_func[mut_info$Func.refGene == "exonic" & mut_info$ExonicFunc.refGene == "nonsynonymous SNV"] <- "Missense"
mut_func[mut_info$Func.refGene == "exonic" & mut_info$ExonicFunc.refGene == "stopgain"] <- "Nonsense"
mut_func[mut_info$Func.refGene == "exonic" & 
           mut_info$ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion")] <- "Frameshift indel"
mut_func[mut_info$Func.refGene == "exonic" & 
           mut_info$ExonicFunc.refGene %in% c("nonframeshift deletion", "nonframeshift insertion")] <- "Inframe indel"

mut_info$mut_func <- mut_func



get_cg_ratio_info <- function(mut_info_selected) {
  

  mut_cg_info <- mut_info_selected %>% 
    group_by(mut_func, cg_type) %>% summarize(count = n()) %>%
    spread(key = cg_type, value = count) %>%
    mutate(CG_ratio = CG / (NonCG + CG))


  CG_ratio_base <- mut_cg_info$CG_ratio[mut_cg_info$mut_func == "Silent"]

  mut_cg_info <- mut_cg_info %>%
    mutate(CG_pV = 1 - pbinom(CG - 1, NonCG + CG, CG_ratio_base)) %>%
    mutate(CG_log_pV = ifelse(CG_pV > 0, -log10(CG_pV), 10))


  gsm_cg_info <- mut_info_selected %>% 
    filter(GSM != "no-change") %>% 
    group_by(GSM, cg_type) %>% summarize(count = n()) %>%
    spread(key = cg_type, value = count) %>%
    mutate(CG_ratio = CG / (NonCG + CG))


  gsm_cg_info <- gsm_cg_info %>%
    mutate(CG_pV = 1 - pbinom(CG - 1, NonCG + CG, CG_ratio_base)) %>%
    mutate(CG_log_pV = ifelse(CG_pV > 0, -log10(CG_pV), 10))


  gsm_cg_info$mut_func <- gsm_cg_info$GSM

  type_order <- c("silent", "missense", "loss-of-function",
                  "exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site", 
                  "intron-retention", "complex")


  cg_info <- rbind(mut_cg_info, gsm_cg_info)


  cg_info$mut_func2 <- factor(cg_info$mut_func,
                              levels = rev(c("Silent", "Missense", "Nonsense", "Inframe indel", "Frameshift indel", "Other",
                                         "exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site", "intron-retention", "complex")),
                              labels = rev(c("Silent", "Missense", "Nonsense", "Inframe indel", "Frameshift indel", "Other",
                                         "Exon skip", "Alternative 5'-ss", "Alternative 3'-ss", "Intron retention", "Complex")))
  cg_info$is_gsm <- ifelse(is.na(cg_info$GSM), "non_gsm", "gsm")


  cg_info_proc <- as.data.frame(cg_info) %>% gather(key = class_statistics, value = value, CG_ratio, CG_log_pV) %>% 
    select(mut_func2, class_statistics, value)

  return(cg_info_proc)
  
}
  

cg_info_proc_master <- data.frame()

cg_info_proc_tmp <- get_cg_ratio_info(mut_info %>% select(GSM, mut_func, cg_type = VogelsteinEtAl_2013))
cg_type_tmp <- rep("Vogelstein (2013)", nrow(cg_info_proc_tmp))
names(cg_type_tmp) <- "Cancer_Gene_Type"
cg_info_proc_master <- rbind(cg_info_proc_master, cbind(cg_info_proc_tmp, cg_type_tmp))

cg_info_proc_tmp <- get_cg_ratio_info(mut_info %>% select(GSM, mut_func, cg_type = YeEtAl_2016))
cg_type_tmp <- rep("Ye (2016)", nrow(cg_info_proc_tmp))
names(cg_type_tmp) <- "Cancer_Gene_Type"
cg_info_proc_master <- rbind(cg_info_proc_master, cbind(cg_info_proc_tmp, cg_type_tmp))

cg_info_proc_tmp <- get_cg_ratio_info(mut_info %>% select(GSM, mut_func, cg_type = LawrenceEtAl_2014))
cg_type_tmp <- rep("Lawrence (2014)", nrow(cg_info_proc_tmp))
names(cg_type_tmp) <- "Cancer_Gene_Type"
cg_info_proc_master <- rbind(cg_info_proc_master, cbind(cg_info_proc_tmp, cg_type_tmp))

cg_info_proc_tmp <- get_cg_ratio_info(mut_info %>% select(GSM, mut_func, cg_type = Cancer_Gene_Census))
cg_type_tmp <- rep("CGC (Feb 2017)", nrow(cg_info_proc_tmp))
names(cg_type_tmp) <- "Cancer_Gene_Type"
cg_info_proc_master <- rbind(cg_info_proc_master, cbind(cg_info_proc_tmp, cg_type_tmp))




ggplot(cg_info_proc_master %>% filter(!(mut_func2 %in% c("Inframe indel", "Frameshift indel", "Other")) & class_statistics == "CG_ratio"),
       aes(x = mut_func2, y = value, fill = mut_func2)) + 
  geom_bar(stat = "identity", position = "dodge")  + 
  coord_flip() +
  labs(x = "", y = "Cancer gene ratio") +
  my_theme() +
  theme(strip.text.x = element_text(size = rel(1.2), angle = 0, hjust = 0),
        panel.spacing.x=unit(1.2, "lines")) +
  facet_grid(cg_type_tmp~., scales = "free") +
  scale_fill_manual(values = splicing_class_colour) +
  # scale_fill_manual(values = c(gsm = "#33a02c", non_gsm = "#b2df8a")) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = FALSE)


ggsave("../figure/cancer_gene_ratio2.tiff", width = 10, height = 10, dpi = 600, units = "cm")


ggplot(cg_info_proc_master %>% filter(!(mut_func2 %in% c("Silent", "Inframe indel", "Frameshift indel", "Other")) & class_statistics == "CG_log_pV"),
       aes(x = mut_func2, y = value, fill = mut_func2)) + 
  geom_bar(stat = "identity", position = "dodge")  + 
  coord_flip() +
  labs(x = "", y = "log10(P-value)") +
  my_theme() +
  theme(strip.text.x = element_text(size = rel(1.2), angle = 0, hjust = 0),
        panel.spacing.x=unit(1.2, "lines")) +
  facet_grid(cg_type_tmp~., scales = "free") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = FALSE)


ggsave("../figure/cancer_gene_pV2.tiff", width = 10, height = 10, dpi = 600, units = "cm")




library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

source("../../../conf/plot_config.R")

cancer_gene <- read.table("../../../db/cancer_gene/cancer_gene.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    select(Gene_Symbol, VogelsteinEtAl_2013) %>%
    filter(VogelsteinEtAl_2013 != "---") 

  

gene_list <- read.table("../../../db/refGene/refGene.txt.gz", sep="\t", stringsAsFactors = FALSE) %>% 
  filter(substring(V2, 0, 3) == "NM_") %>%
  select(Gene_Symbol = V13) %>% distinct(Gene_Symbol) %>% left_join(cancer_gene, key = "Gene_Symbol") 

gene_list$VogelsteinEtAl_2013[is.na(gene_list$VogelsteinEtAl_2013)] <- "NonCG"

print(gene_list)

mut_info <- read.table("../temporary/TCGA.mutation.merged.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  filter(!str_detect(Gene.refGene, "-") & !str_detect(Gene.refGene, ",")) %>%
  inner_join(gene_list, by = c("Gene.refGene" = "Gene_Symbol")) 

print("OK")


mut_func <- rep("Other", nrow(mut_info))
mut_func[mut_info$Func.refGene == "exonic" & mut_info$ExonicFunc.refGene == "synonymous SNV"] <- "Silent"
mut_func[mut_info$Func.refGene == "exonic" & mut_info$ExonicFunc.refGene == "nonsynonymous SNV"] <- "Missense"
mut_func[mut_info$Func.refGene == "exonic" & mut_info$ExonicFunc.refGene == "stopgain"] <- "Nonsense"
mut_func[mut_info$Func.refGene == "exonic" & 
           mut_info$ExonicFunc.refGene %in% c("frameshift deletion", "frameshift insertion")] <- "Frameshift indel"
mut_func[mut_info$Func.refGene == "exonic" & 
           mut_info$ExonicFunc.refGene %in% c("nonframeshift deletion", "nonframeshift insertion")] <- "Inframe indel"

mut_info$mut_func <- mut_func

mut_cg_info <- mut_info %>% 
  group_by(mut_func, VogelsteinEtAl_2013) %>% summarize(count = n()) %>%
  spread(key = VogelsteinEtAl_2013, value = count) %>%
  mutate(Oncogene_ratio = Oncogene / (NonCG + Oncogene + TSG),
         TSG_ratio = TSG / (NonCG + Oncogene + TSG))

TSG_ratio_base <- mut_cg_info$TSG_ratio[mut_cg_info$mut_func == "Silent"]
Oncogene_ratio_base <- mut_cg_info$Oncogene_ratio[mut_cg_info$mut_func == "Silent"]

mut_cg_info <- mut_cg_info %>%
  mutate(Oncogene_pV = 1 - pbinom(Oncogene - 1, NonCG + Oncogene + TSG, Oncogene_ratio_base),
         TSG_pV = 1 - pbinom(TSG - 1, NonCG + Oncogene + TSG, TSG_ratio_base)) %>%
  mutate(Oncogene_log_pV = ifelse(Oncogene_pV > 0, -log10(Oncogene_pV), 10),
         TSG_log_pV = ifelse(TSG_pV > 0, -log10(TSG_pV), 10))


gsm_cg_info <- mut_info %>% 
  filter(GSM != "no-change") %>% 
  group_by(GSM, VogelsteinEtAl_2013) %>% summarize(count = n()) %>%
  spread(key = VogelsteinEtAl_2013, value = count) %>%
  mutate(Oncogene_ratio = Oncogene / (NonCG + Oncogene + TSG),
         TSG_ratio = TSG / (NonCG + Oncogene + TSG))

gsm_cg_info <- gsm_cg_info %>%
  mutate(Oncogene_pV = 1 - pbinom(Oncogene - 1, NonCG + Oncogene + TSG, Oncogene_ratio_base),
         TSG_pV = 1 - pbinom(TSG - 1, NonCG + Oncogene + TSG, TSG_ratio_base)) %>%
  mutate(Oncogene_log_pV = ifelse(Oncogene_pV > 0, -log10(Oncogene_pV), 10),
         TSG_log_pV = ifelse(TSG_pV > 0, -log10(TSG_pV), 10))

gsm_cg_info$mut_func <- gsm_cg_info$GSM

type_order <- c("silent", "missense", "loss-of-function",
                "exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site", 
                "intron-retention", "complex")


cg_info <- rbind(mut_cg_info, gsm_cg_info)

# cg_info$mut_func[cg_info$mut_func == "alternative-5'-splice-site"] <- "alternative-splice-site"
# cg_info$mut_func[cg_info$mut_func == "alternative-3'-splice-site"] <- "alternative-splice-site"


cg_info$mut_func2 <- factor(cg_info$mut_func,
                            levels = rev(c("Silent", "Missense", "Nonsense", "Inframe indel", "Frameshift indel", "Other",
                                       "exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site", "intron-retention", "complex")),
                            labels = rev(c("Silent", "Missense", "Nonsense", "Inframe indel", "Frameshift indel", "Other",
                                       "Exon skip", "Alternative 5'-ss", "Alternative 3'-ss", "Intron retention", "Complex")))
cg_info$is_gsm <- ifelse(is.na(cg_info$GSM), "non_gsm", "gsm")


cg_info_proc <- cg_info %>% gather(key = class_statistics, value = value, Oncogene_ratio, TSG_ratio, Oncogene_log_pV, TSG_log_pV) %>% 
  select(mut_func2, class_statistics, value)

cg_info_proc$VogelsteinEtAl_2013 <- ifelse(cg_info_proc$class_statistics %in% c("Oncogene_ratio", "Oncogene_log_pV"), "Oncogene", "TSG")
cg_info_proc$statistics <- ifelse(cg_info_proc$class_statistics %in% c("Oncogene_ratio", "TSG_ratio"), "ratio", "log_pV")


cg_info_proc$is_gsm <- ifelse(cg_info$mut_func2 %in% 
                                c("Exon skip", "Alternative splice site",
                                  "Intron retention", "Complex"), "gsm", "non_gsm")


ggplot(cg_info_proc %>% filter(!(mut_func2 %in% c("Inframe indel", "Frameshift indel", "Other")) & statistics == "ratio"),
       aes(x = mut_func2, y = value, fill = mut_func2)) + 
  geom_bar(stat = "identity", position = "dodge")  + 
  coord_flip() +
  labs(x = "", y = "Cancer gene ratio") +
  my_theme() +
  theme(strip.text.x = element_text(size = rel(1.2), angle = 0, hjust = 0),
        panel.spacing.x=unit(1.2, "lines")) +
  facet_grid(.~VogelsteinEtAl_2013, scales = "free") +
  scale_fill_manual(values = splicing_class_colour) +
  # scale_fill_manual(values = c(gsm = "#33a02c", non_gsm = "#b2df8a")) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = FALSE)


ggsave("../figure/cancer_gene_ratio.pdf", width = 8, height = 2.2)


ggplot(cg_info_proc %>% filter(!(mut_func2 %in% c("Silent", "Inframe indel", "Frameshift indel", "Other")) & statistics == "log_pV"),
       aes(x = mut_func2, y = value, fill = mut_func2)) + 
  geom_bar(stat = "identity", position = "dodge")  + 
  coord_flip() +
  labs(x = "", y = "log10(P-value)") +
  my_theme() +
  theme(strip.text.x = element_text(size = rel(1.2), angle = 0, hjust = 0),
        panel.spacing.x=unit(1.2, "lines")) +
  facet_grid(.~VogelsteinEtAl_2013, scales = "free") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = FALSE)


ggsave("../figure/cancer_gene_pV.pdf", width = 8, height = 2.0)


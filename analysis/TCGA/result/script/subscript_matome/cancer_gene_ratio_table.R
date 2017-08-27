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

cg_count_table <- data.frame(cg_type = c(), cg_num = c(), non_cg_num = c(), cg_ratio = c())


get_cg_stat <- function(mut_info) {

    # Number of CG SAVS
    a <- mut_info %>% filter(GSM != "no-change" & cg_type == "CG")
    n1 <- nrow(a)

    # Number of CG SAVs disrupting non-canonical splice sites
    a <- mut_info %>% filter(GSM != "no-change" & cg_type == "CG" & grepl("disruption", Mutation_Type) & Is_Canonical == "non-canonical")
    n2 <- nrow(a) 

    # Number of CG SAVs creationg novel splice sites
    a <- mut_info %>% filter(GSM != "no-change" & cg_type == "CG" & grepl("creation", Mutation_Type))
    n3 <- nrow(a)

    # Number of samples with at least one CG SAVs
    a <- mut_info %>% filter(GSM != "no-change" & cg_type == "CG") %>% select(Sample_Name) %>% distinct()
    n4 <- nrow(a)

    # Number of samples with at least one CG SAVs other than those disrupting canonical splices sites
    a <- mut_info %>% filter(GSM != "no-change" & cg_type == "CG" & (grepl("creation", Mutation_Type) | Is_Canonical == "non-canonical")) %>%
           select(Sample_Name) %>% distinct()
    n5 <- nrow(a)

    # Ratio of SAVs to total LOF variants
    b1 <- mut_info %>% filter((GSM != "no-change" | ExonicFunc.refGene %in% c("stopgain", "frameshift deletion", "frameshift insertion")) & cg_type == "CG")
    b2 <- mut_info %>% filter(GSM != "no-change" & cg_type == "CG")
    n6 <- nrow(b2) / nrow(b1)

    return(c(n1, n2, n3, n4, n5, n6))
}

V <- get_cg_stat(mut_info %>% mutate(cg_type = VogelsteinEtAl_2013))
V <- rbind(V, get_cg_stat(mut_info %>% mutate(cg_type = YeEtAl_2016)))
V <- rbind(V, get_cg_stat(mut_info %>% mutate(cg_type = LawrenceEtAl_2014)))
V <- rbind(V, get_cg_stat(mut_info %>% mutate(cg_type = Cancer_Gene_Census))) 

print(V)
V <- cbind(c("VogelsteinEtAl_2013", "YeEtAl_2016", "LawrenceEtAl_2014", "Cancer_Gene_Census"), V)

colnames(V) <- c("Cancer_Gene_Type", "CG_SAV_COUNT", "CG_SAV_ND_COUNT", "CG_SAV_C_COUNT", "CG_SAV_SAMPLE_COUNT", "CG_SAV_C_ND_SAMPLE_COUNT", "CG_SAV_RATIO_LOF")

write.table(V, "../table/cancer_gene_ratio_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

 



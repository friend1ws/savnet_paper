library(dplyr)
library(ggplot2)
library(VennDiagram)


target_sample <- substring(read.table("../data/target_sample_list_JungEtAl.txt", 
                            header = FALSE, stringsAsFactors = FALSE)$V1, 1, 12)

cancer_gene <- read.table("/home/yshira/mysoftware/sv_utils/cancer_gene/cancer_gene.txt", sep = "\t", header = FALSE)
cancer_gene_list <- as.character(cancer_gene[cancer_gene$V5 == "CG", "V1"])


JungEtAl_RS <- read.table("../data/ng.3414-S2.ST3.proc.txt", sep = "\t", 
                          header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  filter(Distance.to.5..donor.splice.site %in% c(-2, 1, -1, -2) |
         Distance.to.3..acceptor.splice.site %in% c(1, -1, -2)) %>%
  filter(Sample %in% target_sample) %>%
  mutate(Comp_Key = paste(Sample, sub("chr", "", Chromosome), Position, sep = "_"))

JungEtAl_RS_CG <- JungEtAl_RS %>% filter(Gene %in% cancer_gene_list)


omega <- read.table("../../matome/omega.motif_summary.txt", sep = "\t",
                    header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  filter(substring(Sample_Name, 1, 12) %in% target_sample) %>%
  filter((Mutation_Type == "splicing donor disruption" & Rel_Pos %in% c(1, 2, 3, 4, 5)) |
          Mutation_Type == "splicing acceptor disruption" & Rel_Pos %in% c(5, 6, 7)) %>%
  mutate(Comp_Key = paste(substring(Sample_Name, 1, 12), 
                          unlist(lapply(strsplit(Mutation_Key, ','), function(x) {paste(x[1], x[2], sep = "_")})),
                          sep = "_"))

omega_CG <- omega %>% filter(Gene_Symbol %in% cancer_gene_list)


All_comp <- list(GSM = omega$Comp_Key, JungEtAl = JungEtAl_RS$Comp_Key)

venn.diagram(All_comp, 
             fill = c("dodgerblue", "goldenrod1"),
             filename = "../output/JungEtAl_comp.png",
             height = 3000,
             width = 3000,
             margin = 0.05,
             main = "all genes",
             main.cex = 1.5)

CG_comp <- list(GSM = omega_CG$Comp_Key, JungEtAl = JungEtAl_RS_CG$Comp_Key)

venn.diagram(CG_comp, 
             fill = c("dodgerblue", "goldenrod1"),
             filename = "../output/JungEtAl_comp_CG.png",
             height = 3000,
             width = 3000,
             margin = 0.05,
             main = "cancer genes",
             main.cex = 1.5)



library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

# source("../../../conf/plot_config.R")
source("subscript_matome/sav_profile_function.R")

# read_num_info <- read.table("../temporary/TCGA.motif_read_num.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
sav_result <- read.table("../../output/rescue/TCGA.savnet.with_rescued.result.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
    filter(Mutation_Type %in% c("splicing donor disruption", "splicing acceptor disruption"))


sav_result$Mutation_Type <- factor(sav_result$Mutation_Type,
                                levels = c("splicing donor disruption", "splicing acceptor disruption"),
                                labels = c("Donor", "Acceptor"))

sav_result$Is_Canonical <- factor(sav_result$Is_Canonical,
                                  levels = c("canonical", "non-canonical"),
                                  labels = c("Canonical", "Non-canonical"))


refGene <- read.table("../../../db/refGene/refGene.txt.gz", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# gene2ref_tmp <- read.table("../db/gene2ref.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene2ref_tmp <- read.table("../../../db/uniprot/gene2ref.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene2ref <- gene2ref_tmp$V2
names(gene2ref) <- gene2ref_tmp$V1




# filter CDKN2A NM_058195
check_CDKN2A <- function(x) {
  motif_pos_end <- unlist(lapply(strsplit(x$Splicing_Key, "-"), '[', 2))
  return(x$Gene_Symbol == "CDKN2A" & motif_pos_end > 21975132)
}

sav_result <- sav_result[!check_CDKN2A(sav_result), ]


mutation_count <- sav_result %>% select(Gene_Symbol, Sample_Name, Mutation_Key, Mutation_Type, Motif_Pos, Is_Canonical, YeEtAl_2016) %>% 
                                        distinct() %>% 
                                        # filter(Gene_Symbol %in% names(gene2ref)) %>%
                                        # filter(YeEtAl_2016 == "CG") %>%
                                        group_by(Mutation_Key, Gene_Symbol, Mutation_Type, Motif_Pos, Is_Canonical) %>% 
                                        summarize(count = n()) %>%
                                        mutate(Chr = unlist(lapply(strsplit(Mutation_Key, split = ","), '[[', 1)),
                                               Pos = unlist(lapply(strsplit(Mutation_Key, split = ","), '[[', 2)),
                                               Ref = unlist(lapply(strsplit(Mutation_Key, split = ","), '[[', 3)),
                                               Alt = unlist(lapply(strsplit(Mutation_Key, split = ","), '[[', 4))) %>%
                                        mutate(Sub = paste(Ref, Alt, sep = ">")) %>%
                                        mutate(Sub_Count = paste(Pos, Sub, Is_Canonical, count, sep = ","))
                                        # mutate(Sub_Count = paste(Sub, count, sep = ";")) 



hotspot_summary <- mutation_count %>% 
  group_by(Gene_Symbol, Motif_Pos) %>%
  summarize(Total_Count = sum(count),
            Substitution_Count = paste(Sub_Count, collapse = ";")) %>% 
  filter(Total_Count >= 5) %>%
  arrange(Gene_Symbol, Motif_Pos)



motif_infos <- lapply(1:nrow(hotspot_summary), 
                     function(i) {
                       get_motif_info(as.character(hotspot_summary$Motif_Pos[i]), 
                                      refGene,  
                                      as.character(gene2ref[as.character(hotspot_summary$Gene_Symbol[i])]))})


hotspot_summary$Motif_Type <- unlist(lapply(motif_infos, '[', 1))
hotspot_summary$Exon_Num <- unlist(lapply(motif_infos, '[', 2))


hotspot_summary$Motif_Type <- factor(hotspot_summary$Motif_Type,
                                     levels = c("donor", "acceptor"),
                                     labels = c("Donor", "Acceptor"))

motif_pos_sp1 <- strsplit(hotspot_summary$Motif_Pos, ":")
motif_pos_sp2 <- strsplit(unlist(lapply(motif_pos_sp1, '[', 2)), ",")
motif_pos_sp3 <- strsplit(unlist(lapply(motif_pos_sp2, '[', 1)), "-")

tchr <- unlist(lapply(motif_pos_sp1, '[', 1))
tstart <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 1)))
tend <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 2)))
tstrand <- unlist(lapply(motif_pos_sp2, '[', 2))

hotspot_summary$Chr <- paste("chr", tchr, sep = "") 
hotspot_summary$Start <- tstart 
hotspot_summary$End <- tend
hotspot_summary$Strand <- tstrand

hotspot_summary$Transcript <- gene2ref[hotspot_summary$Gene_Symbol]

hotspot_summary <- as.data.frame(hotspot_summary) %>% select(Gene_Symbol, Transcript, Chr, Start, End, Strand, Motif_Type, Exon_Num, Total_Count, Substitution_Count)

colnames(hotspot_summary) <- c("Gene", "Transcript", "Chr", "Start", "End", "Strand", "Motif", "Exon number", "Total count", "Each substitution (Pos, Substitution, Canonical or not, Count)")

write.table(hotspot_summary, "../table/TableS7.txt", quote = FALSE, row.names = FALSE, sep = "\t")




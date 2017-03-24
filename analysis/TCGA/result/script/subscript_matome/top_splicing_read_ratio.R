library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

source("../../../conf/plot_config.R")

read_num_info <- read.table("../temporary/TCGA.motif_read_num.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)
refGene <- read.table("../../../db/refGene/refGene.txt.gz", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene2ref_tmp <- read.table("../db/gene2ref.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene2ref <- gene2ref_tmp$V2
names(gene2ref) <- gene2ref_tmp$V1


# get intron pos for opposite intron retention
get_motif_info <- function(motif_pos, ref_gene_id) {
  
  motif_pos_sp1 <- strsplit(motif_pos, ":")[[1]]
  motif_pos_sp2 <- strsplit(motif_pos_sp1[2], ",")[[1]]
  motif_pos_sp3 <- strsplit(motif_pos_sp2[1], "-")[[1]]
  
  tchr <- motif_pos_sp1[1]
  tstart <- as.numeric(motif_pos_sp3[1])
  tend <- as.numeric(motif_pos_sp3[2])
  tdir <- motif_pos_sp2[2]
  
  target_gene_info <- refGene %>% filter(V2 == ref_gene_id)
  exon_starts <- as.numeric(strsplit(as.character(target_gene_info[10]), split=",")[[1]]) + 1
  exon_ends <- as.numeric(strsplit(as.character(target_gene_info[11]), split=",")[[1]])
  
  exon_ind1 <- which(exon_starts >= tstart & exon_starts <= tend)
  exon_ind2 <- which(exon_ends >= tstart & exon_ends <= tend)  
  
  if (length(exon_ind1) > 0) {
    c(ifelse(tdir == "+", "acceptor", "donor"), 
      ifelse(tdir == "+", exon_ind1, length(exon_starts) - exon_ind1 + 1))
  } else if (length(exon_ind2) > 0) {
    c(ifelse(tdir == "+", "donor", "acceptor"), 
      ifelse(tdir == "+", exon_ind2, length(exon_starts) - exon_ind2 + 1))  
  } else {
    return(c(NA, NA))
  }
  
}



motif2count <- read_num_info %>% 
  select(Sample_Name, Gene_Symbol, Motif_Pos) %>% 
  distinct() %>% 
  group_by(Motif_Pos, Gene_Symbol) %>% 
  summarize(count = n()) %>% 
  arrange(desc(count)) %>% 
  filter(count >= 8)


Ds <- data.frame()
for (i in 1:nrow(motif2count)) {
  
  read_num_info_filt <- read_num_info %>% 
    filter(Motif_Pos == as.character(motif2count[i, 1])) %>%
    mutate(Splicing_Key2 = paste(Splicing_Key, Splicing_Class, sep = " ")) %>%
    mutate(n_Supporting_Read_Num = Supporting_Read_Num / Weight) %>% 
    select(Sample_Name, Splicing_Key2, Mutation_Key, n_Supporting_Read_Num) %>% 
    spread(key = Splicing_Key2, value = n_Supporting_Read_Num)
  
  mut_count <- table(read_num_info_filt$Mutation_Key)
  mut_id <- factor(read_num_info_filt$Mutation_Key,
                   levels = names(mut_count)[order(mut_count, decreasing = TRUE)],
                   labels = 1:length(mut_count))
  
  N <- nrow(read_num_info_filt)
  
  int_id_1 <- grep(" intron-retention", colnames(read_num_info_filt))
  int_id_2 <- grep(" opposite-side-intron-retention", colnames(read_num_info_filt))
  
  if (length(int_id_1) >= 1 & length(int_id_2) >= 1) {
    read_num_info_filt[, int_id_1] <- pmax(read_num_info_filt[,int_id_1], read_num_info_filt[, int_id_2])
    read_num_info_filt[, int_id_2] <- 0
  }
  
  vec_mat <- read_num_info_filt[,3:ncol(read_num_info_filt), drop = FALSE]
  mean_count_vec <- apply(vec_mat, 2, function(x) {mean(x, trim = 0.1)})
  order_vec <- order(mean_count_vec, decreasing = TRUE)
  
  top_ind <- order_vec[1]
  top_splice <- names(mean_count_vec)[top_ind]
  first_vec <- vec_mat[, top_ind]
  
  sec_ind <- 0
  sec_splice <- "---"
  if (length(order_vec) >= 2 & mean_count_vec[order_vec[2]] > 0) {
    vec_mat[,top_ind] <- 0 
    second_vec <- as.numeric(apply(vec_mat, 1, max))
    sec_splice <- names(mean_count_vec)[sec_ind]
    rel_vector <- first_vec / (first_vec + second_vec)
    rel_vector[is.na(rel_vector)] <- 0.5
    # rel_vector <- pmin(rel_vector, 2)
  } else {
    rel_vector <- rep(1, N)
  }
  
  Ds <- rbind(Ds, data.frame(Motif_Pos = rep(as.character(motif2count[i, 1]), N),
                             Gene_Symbol = rep(as.character(motif2count[i, 2]), N),
                             Mut_Num = rep(N, N),
                             Top_Splice = rep(strsplit(top_splice, ' ')[[1]][2],  N),
                             Rel_Count = rel_vector,
                             Mut_ID = mut_id))
  
}


Ds$Top_Splice[Ds$Top_Splice == "intronic-alternative-3'-splice-site"] <- "alternative-3'-splice-site"
Ds$Top_Splice[Ds$Top_Splice == "intronic-alternative-5'-splice-site"] <- "alternative-5'-splice-site"

Ds$Top_Splice <- factor(Ds$Top_Splice,
                        levels = c("exon-skip", "alternative-5'-splice-site", 
                                   "alternative-3'-splice-site", "intron-retention"),
                        labels = c("Exon skip", "Alternative 5'-ss", "Alternative 3'-ss", "Intron retention"))

motif2mut_splice_info <- Ds %>% select(Motif_Pos, Mut_Num, Top_Splice) %>% distinct()




DDs <- Ds %>% select(Motif_Pos, Gene_Symbol) %>% distinct()
 
 
motif_infos <- lapply(1:nrow(DDs), 
                     function(i) {
                       get_motif_info(as.character(DDs$Motif_Pos[i]), 
                                      as.character(gene2ref[as.character(DDs$Gene_Symbol[i])]))})



motif2info <- data.frame(Motif_Pos = DDs$Motif_Pos, 
                         Gene_Symbol = factor(DDs$Gene_Symbol, 
                                              levels = c("TP53", "PIK3R1", "CDKN2A", "GATA3", "MET", "MIEN1")),
                         Motif_Type = factor(unlist(lapply(motif_infos, '[', 1)),
                                             levels = c("acceptor", "donor")),
                         Exon_Num = as.numeric(unlist(lapply(motif_infos, '[', 2)))) %>%
  arrange(Gene_Symbol, Exon_Num, Motif_Type)

motif2info <- motif2info %>% full_join(motif2mut_splice_info, by = "Motif_Pos")




# motif_table <- data.frame(x = c(0.2, 0.5, 0.9, 1.1, 1.8), 
#                           y = c(0, 0, 0, 0, 0), 
#                           label = c("Gene", "Motif", "Exon", "Major splicing", "#SAVs"),
#                           hjust = c(0, 0, 1, 0, 1))


# current_y <- 1
# for (i in 1:nrow(motif2info)) {
  
#   motif_table <- rbind(motif_table, 
#         data.frame(x = c(0.2, 0.5, 0.9, 1.1, 1.8),
#                    y = rep(current_y, 5),
#                    label = c(as.character(motif2info[i, "Gene_Symbol"]),
#                              as.character(motif2info[i, "Motif_Type"]),
#                              as.character(motif2info[i, "Exon_Num"]),
#                              as.character(motif2info[i, "Top_Splice"]),
#                              as.character(motif2info[i, "Mut_Num"])),
#                    hjust = c(0, 0, 1, 0, 1)))
  
#   current_y <- current_y + 1
  
# }


  
# p1 <- ggplot(motif_table %>% filter(y != 0), aes(x = x, y = y, label = label, hjust = hjust)) + 
#   geom_text(size = 3) + 
#   scale_y_continuous(trans = "reverse", breaks = unique(motif_table$y)) +
#   theme_bw() +
#   labs(x = "", y = "") + 
#   theme(panel.grid = element_blank(),
#         panel.border = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())




        

DDDs <- Ds %>% left_join(motif2mut_splice_info, by = "Motif_Pos") %>%
  left_join(motif2info, by = "Motif_Pos") %>% 
  arrange(Gene_Symbol.x, Exon_Num, Top_Splice.x)




key <- paste(DDDs$Gene_Symbol.x, " exon ", 
             DDDs$Exon_Num, ", ", 
             DDDs$Motif_Type, " (", 
             DDDs$Top_Splice.x, ")" , sep = "")

# paste(Ds$Gene_Symbol, Ds$Motif_Pos, Ds$Top_Splice, Ds$Mut_Num)
DDDs$Key <- factor(key, levels = rev(unique(key)))




# ggsave("../figure/top_splicing_read_ratio.pdf", width = 8, height = 8)
ggplot(DDDs %>% filter(Mut_ID %in% as.character(1:12)), aes(x = Key, y = Rel_Count, colour = Mut_ID)) + 
  geom_jitter(width = 0.15, height = 0.02, size = 0.8) + 
  coord_flip() +
  my_theme() +
  # theme_minimal() +
  # theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank()) +
  theme(axis.ticks = element_blank(),
        panel.grid.major.y = element_line(linetype = "longdash", colour = "grey60")) +
  labs(x = "", y = "Fraction of the most frequent type of abnormal splicing") +
  # labs(x = "", y = "") +
    guides(colour = FALSE) 


ggsave("../figure/top_splicing_read_ratio.tiff", width = 18, height = 9, dpi = 600, units = "cm")





##########
# example

motif2first_second_read_count <- function(motif_pos_str) {
  

  read_num_info_filt <- read_num_info %>% 
    filter(Motif_Pos == motif_pos_str) %>%
    mutate(Splicing_Key2 = paste(Splicing_Class, " (", Splicing_Key, ")", sep = "")) %>%
    # mutate(n_Supporting_Read_Num = Supporting_Read_Num / Weight) %>% 
    select(Sample_Name, Splicing_Key2, Mutation_Key, Supporting_Read_Num) %>% 
    spread(key = Splicing_Key2, value = Supporting_Read_Num)

  mut_count <- table(read_num_info_filt$Mutation_Key)
  mut_id <- factor(read_num_info_filt$Mutation_Key,
                   levels = names(mut_count)[order(mut_count, decreasing = TRUE)],
                   labels = 1:length(mut_count))

  N <- nrow(read_num_info_filt)

  int_id_1 <- grep("intron-retention", colnames(read_num_info_filt))
  int_id_2 <- grep("opposite-side-intron-retention", colnames(read_num_info_filt))
  int_id_1 <- setdiff(int_id_1, int_id_2)
  
  if (length(int_id_1) >= 1 & length(int_id_2) >= 1) {
    read_num_info_filt[, int_id_1] <- pmax(read_num_info_filt[,int_id_1], read_num_info_filt[, int_id_2])
    read_num_info_filt[, int_id_2] <- 0
  }

  vec_mat <- read_num_info_filt[,3:ncol(read_num_info_filt), drop = FALSE]
  mean_count_vec <- apply(vec_mat, 2, function(x) {mean(x, trim = 0.1)})
  order_vec <- order(mean_count_vec, decreasing = TRUE)

  top_ind <- order_vec[1]
  sec_ind <- order_vec[2]
  top_splice <- names(mean_count_vec)[top_ind]
  sec_splice <- names(mean_count_vec)[sec_ind]
  
  top_splice <- sub("exon-skip", "Exon skip", top_splice)
  top_splice <- sub("intronic-", "", top_splice)
  top_splice <- sub("opposite-side-", "", top_splice)
  top_splice <- sub("intron-retention", "Intron retention", top_splice)
  top_splice <- sub("alternative-5'-splice-site", "Alternative 5'-ss", top_splice)
  top_splice <- sub("alternative-3'-splice-site", "Alternative 3'-ss", top_splice)
  
  sec_splice <- sub("exon-skip", "Exon skip", sec_splice)
  sec_splice <- sub("intronic-", "", sec_splice)
  sec_splice <- sub("opposite-side-", "", sec_splice)
  sec_splice <- sub("intron-retention", "Intron retention", sec_splice)
  sec_splice <- sub("alternative-5'-splice-site", "Alternative 5'-ss", sec_splice)
  sec_splice <- sub("alternative-3'-splice-site", "Alternative 3'-ss", sec_splice)
  
  first_vec <- vec_mat[, top_ind]
  second_vec <- vec_mat[, sec_ind]
  
  
  temp_for_gene_symbol <- read_num_info %>% filter(Motif_Pos == motif_pos_str)
  gene_symbol <- unique(temp_for_gene_symbol$Gene_Symbol)
  motif_info <- get_motif_info(motif_pos_str, as.character(gene2ref[gene_symbol])) 
  plot_title <- paste(gene_symbol, " exon ", motif_info[2], ", ", motif_info[1], sep = "")
  
  df <- data.frame(mut_id = mut_id, first_vec = first_vec, second_vec = second_vec)
  
  p <- ggplot(df, aes(x = first_vec, y = second_vec, colour = mut_id)) +
    geom_point(size = 1) +
    my_theme() +
    ggtitle(plot_title) +
    labs(x = top_splice, y = sec_splice) +
    guides(colour = FALSE)
  
  return(p)

}


p_TP53 <- motif2first_second_read_count("17:7579306-7579314,-")
p_CDKN2A <- motif2first_second_read_count("9:21971207-21971213,-")
p_GATA3 <- motif2first_second_read_count("10:8111430-8111436,+")


plot_grid(p_TP53, p_CDKN2A, p_GATA3, ncol = 3)



ggsave("../figure/top_splicing_read_ratio_example.tiff", width = 20, height = 7, dpi = 600, units = "cm")




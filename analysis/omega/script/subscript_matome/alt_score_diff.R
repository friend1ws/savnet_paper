library(dplyr)
library(spliceSites)
library(stringr)
library(ggplot2)

mes <- load.maxEnt()
hb <- load.hbond()

get_pos_df <- function(motif_pos) {
  
  motif_pos_sp1 <- strsplit(motif_pos, ":")
  tchr <- paste("chr", unlist(lapply(motif_pos_sp1, '[', 1)), sep = "")
  
  motif_pos_sp2 <- strsplit(unlist(lapply(motif_pos_sp1, '[', 2)), ",")
  tstrand <- unlist(lapply(motif_pos_sp2, '[', 2))
  
  motif_pos_sp3 <- strsplit(unlist(lapply(motif_pos_sp2, '[', 1)), "-")
  tstart <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 1)))
  tend <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 2)))
  
  return(data.frame(chr = tchr, start = tstart, end = tend, strand = tstrand))
  
}



##########


pos_info <- read.table("../matome/omega.alt_junc.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)

mes_info <- read.table("../matome/omega.motif_summary.mes.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)

mut_info <- inner_join(pos_info, mes_info, key = c("Cancer_Type", "Sample_Name", "Mutation_Key"))

mut_info$Is_SNV <- unlist(lapply(strsplit(mut_info$Mutation_Key, ','), function(x) {str_length(x[3]) == 1 & str_length(x[4]) == 1}))

mut_info <- mut_info %>% filter(Is_SNV == TRUE)

##########
# donor creation

mut_info_dc <- mut_info %>% filter(Mutation_Type == "splicing donor creation")
t_chr_dc <- unlist(lapply(strsplit(mut_info_dc$Mutation_Key, ','), function(x) {paste("chr", x[1], sep = "")}))
t_strand_dc <- unlist(lapply(strsplit(mut_info_dc$Motif_Pos, ','), function(x) {x[2]}))
  
t_start_dc <- rep(0, nrow(mut_info_dc))
t_end_dc <- rep(0, nrow(mut_info_dc))

t_start_dc[mut_info_dc$Exon_Strand == "+"] <-
  mut_info_dc$Exon_End[mut_info_dc$Exon_Strand == "+"] - 2
t_end_dc[mut_info_dc$Exon_Strand == "+"] <-
  mut_info_dc$Exon_End[mut_info_dc$Exon_Strand == "+"] + 8

t_start_dc[mut_info_dc$Exon_Strand == "-"] <-
  mut_info_dc$Exon_Start[mut_info_dc$Exon_Strand == "-"] - 7
t_end_dc[mut_info_dc$Exon_Strand == "-"] <-
  mut_info_dc$Exon_Start[mut_info_dc$Exon_Strand == "-"] + 3


gr_dc <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = t_chr_dc, start = t_start_dc, end = t_end_dc, strand = t_strand_dc), ignore.strand = FALSE)

seq_org_dc <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr_dc)

mes_org_dc <- unlist(lapply(as.character(seq_org_dc), function(x) {score5(mes, x, pos = 3)}))

hb_score_org_dc <- hbond(hb, as.character(seq_org_dc), 3)


##########
# acceptor creation

mut_info_ac <- mut_info %>% filter(Mutation_Type == "splicing acceptor creation")
t_chr_ac <- unlist(lapply(strsplit(mut_info_ac$Mutation_Key, ','), function(x) {paste("chr", x[1], sep = "")}))
t_strand_ac <- unlist(lapply(strsplit(mut_info_ac$Motif_Pos, ','), function(x) {x[2]}))

t_start_ac <- rep(0, nrow(mut_info_ac))
t_end_ac <- rep(0, nrow(mut_info_ac))

t_start_ac[mut_info_ac$Exon_Strand == "+"] <- mut_info_ac$Exon_Start[mut_info_ac$Exon_Strand == "+"] - 19
t_end_ac[mut_info_ac$Exon_Strand == "+"] <- mut_info_ac$Exon_Start[mut_info_ac$Exon_Strand == "+"] + 3

t_start_ac[mut_info_ac$Exon_Strand == "-"] <- mut_info_ac$Exon_End[mut_info_ac$Exon_Strand == "-"] - 2
t_end_ac[mut_info_ac$Exon_Strand == "-"] <- mut_info_ac$Exon_End[mut_info_ac$Exon_Strand == "-"] + 20


gr_ac <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = t_chr_ac, start = t_start_ac, end = t_end_ac, strand = t_strand_ac), ignore.strand = FALSE)

seq_org_ac <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr_ac)

mes_org_ac <- unlist(lapply(as.character(seq_org_ac), function(x) {score3(mes, x, pos = 20)}))




score_creation <- data.frame(diff_score = c(mut_info_dc$mes_wt - mes_org_dc, 
                                        mut_info_dc$mes_mt - mes_org_dc, 
                                        mut_info_dc$hbond_wt - hb_score_org_dc, 
                                        mut_info_dc$hbond_mt - hb_score_org_dc, 
                                        mut_info_ac$mes_wt - mes_org_ac,
                                        mut_info_ac$mes_mt - mes_org_ac),
                         is_mutation = factor(c(rep("wt", nrow(mut_info_dc)),
                                             rep("mt", nrow(mut_info_dc)),
                                             rep("wt", nrow(mut_info_dc)),
                                             rep("mt", nrow(mut_info_dc)),
                                             rep("wt", nrow(mut_info_ac)),
                                             rep("mt", nrow(mut_info_ac))),
                                           levels = c("wt", "mt")),
                         score_type = factor(c(rep("mes", 2 * nrow(mut_info_dc)),
                                               rep("hbond", 2 * nrow(mut_info_dc)),
                                               rep("mes", 2 * nrow(mut_info_ac))),
                                             levels = c("mes", "hbond")),
                         is_dc = factor(c(rep("donor", 4 * nrow(mut_info_dc)),
                                          rep("acceptor", 2 * nrow(mut_info_ac))),
                                        levels = c("donor", "acceptor"))
                         )


ggplot(score_creation %>% filter(score_type == "mes"), aes(x = is_mutation, y = diff_score, fill = is_mutation)) + 
  geom_boxplot() +
  ylim(c(-25, 25)) +
  labs(x = "", y = "Diff. of MaxEnt score between \nalternative and authentic splice sites") +
  theme_bw() +
  facet_grid(.~is_dc) +
  scale_fill_brewer(palette = "Accent") +
  guides(fill = FALSE) 


ggsave("../matome/diff_mes_creation.png", width = 6, height = 4)


ggplot(score_creation %>% filter(score_type == "hbond"), aes(x = is_mutation, y = diff_score, fill = is_mutation)) + 
  geom_boxplot() +
  ylim(c(-20, 20)) +
  labs(x = "", y = "Diff. of H-bond score between \nalternative and authentic splice sites") +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  guides(fill = FALSE) 

ggsave("../matome/diff_hb_creation.png", width = 3, height = 4)

  
                                    
########
# donor disruption

mut_info_dd <- mut_info %>% filter(Mutation_Type == "splicing donor disruption")
t_chr_dd <- unlist(lapply(strsplit(mut_info_dd$Mutation_Key, ','), function(x) {paste("chr", x[1], sep = "")}))
t_strand_dd <- unlist(lapply(strsplit(mut_info_dd$Motif_Pos, ','), function(x) {x[2]}))

t_start_dd <- rep(0, nrow(mut_info_dd))
t_end_dd <- rep(0, nrow(mut_info_dd))

t_start_dd[mut_info_dd$Exon_Strand == "+"] <- mut_info_dd$Junc_Pos[mut_info_dd$Exon_Strand == "+"] - 3
t_end_dd[mut_info_dd$Exon_Strand == "+"] <- mut_info_dd$Junc_Pos[mut_info_dd$Exon_Strand == "+"] + 7

t_start_dd[mut_info_dd$Exon_Strand == "-"] <- mut_info_dd$Junc_Pos[mut_info_dd$Exon_Strand == "-"] - 7
t_end_dd[mut_info_dd$Exon_Strand == "-"] <- mut_info_dd$Junc_Pos[mut_info_dd$Exon_Strand == "-"] + 3


gr_dd <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = t_chr_dd, start = t_start_dd, end = t_end_dd, strand = t_strand_dd), ignore.strand = FALSE)

seq_alt_dd <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr_dd)

mes_alt_dd <- unlist(lapply(as.character(seq_alt_dd), function(x) {score5(mes, x, pos = 3)}))

hb_score_alt_dd <- hbond(hb, as.character(seq_alt_dd), 3)

##########
# acceptor disruption

mut_info_ad <- mut_info %>% filter(Mutation_Type == "splicing acceptor disruption")
t_chr_ad <- unlist(lapply(strsplit(mut_info_ad$Mutation_Key, ','), function(x) {paste("chr", x[1], sep = "")}))
t_strand_ad <- unlist(lapply(strsplit(mut_info_ad$Motif_Pos, ','), function(x) {x[2]}))

t_start_ad <- rep(0, nrow(mut_info_ad))
t_end_ad <- rep(0, nrow(mut_info_ad))

t_start_ad[mut_info_ad$Exon_Strand == "+"] <- mut_info_ad$Junc_Pos[mut_info_ad$Exon_Strand == "+"] - 19
t_end_ad[mut_info_ad$Exon_Strand == "+"] <- mut_info_ad$Junc_Pos[mut_info_ad$Exon_Strand == "+"] + 3

t_start_ad[mut_info_ad$Exon_Strand == "-"] <- mut_info_ad$Junc_Pos[mut_info_ad$Exon_Strand == "-"] - 3
t_end_ad[mut_info_ad$Exon_Strand == "-"] <- mut_info_ad$Junc_Pos[mut_info_ad$Exon_Strand == "-"] + 19


gr_ad <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = t_chr_ad, start = t_start_ad, end = t_end_ad, strand = t_strand_ad), ignore.strand = FALSE)

seq_alt_ad <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr_ad)

mes_alt_ad <- unlist(lapply(as.character(seq_alt_ad), function(x) {score3(mes, x, pos = 20)}))



score_disruption <- data.frame(diff_score = c(mes_alt_dd - mut_info_dd$mes_wt, 
                                              mes_alt_dd - mut_info_dd$mes_mt,
                                              hb_score_alt_dd - mut_info_dd$hbond_wt,
                                              hb_score_alt_dd - mut_info_dd$hbond_mt,
                                              mes_alt_ad - mut_info_ad$mes_wt,
                                              mes_alt_ad - mut_info_ad$mes_mt),
                             is_mutation = factor(c(rep("wt", nrow(mut_info_dd)),
                                                    rep("mt", nrow(mut_info_dd)),
                                                    rep("wt", nrow(mut_info_dd)),
                                                    rep("mt", nrow(mut_info_dd)),
                                                    rep("wt", nrow(mut_info_ad)),
                                                    rep("mt", nrow(mut_info_ad))),
                                                  levels = c("wt", "mt")),
                             score_type = factor(c(rep("mes", 2 * nrow(mut_info_dd)),
                                                   rep("hbond", 2 * nrow(mut_info_dd)),
                                                   rep("mes", 2 * nrow(mut_info_ad))),
                                                 levels = c("mes", "hbond")),
                             is_dc = factor(c(rep("donor", 4 * nrow(mut_info_dd)),
                                              rep("acceptor", 2 * nrow(mut_info_ad))),
                                            levels = c("donor", "acceptor"))
)


ggplot(score_disruption %>% filter(score_type == "mes"), 
       aes(x = is_mutation, y = diff_score, fill = is_mutation)) + 
  geom_boxplot() +
  ylim(c(-25, 25)) +
  labs(x = "", y = "Diff. of MaxEnt score between \nalternative and authentic splice sites") +
  theme_bw() +
  facet_grid(.~is_dc) +
  scale_fill_brewer(palette = "Accent") +
  guides(fill = FALSE) 


ggsave("../matome/diff_mes_disruption.png", width = 6, height = 4)


ggplot(score_disruption %>% filter(score_type == "hbond"), aes(x = is_mutation, y = diff_score, fill = is_mutation)) + 
  geom_boxplot() +
  ylim(c(-20, 20)) +
  labs(x = "", y = "Diff. of H-bond score between \nalternative and authentic splice sites") +
  theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  guides(fill = FALSE) 

ggsave("../matome/diff_hb_disruption.png", width = 3, height = 4)


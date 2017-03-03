library(dplyr)
library(spliceSites)
library(ggplot2)
library(cowplot)

source("../../../conf/plot_config.R")

mes<-load.maxEnt()
hb<-load.hbond()


# exon_size_d <- 3
# intron_size_d <- 6
# intron_size_a <- 6
# exon_size_a <- 1

snv_info <- read.table("../temporary/TCGA.savnet.allele_count.summary.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "") %>%
  filter(Ref_Mut != "-" & Alt_Mut != "-")

snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "A"] <- "S"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "T"] <- "A"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "S"] <- "T"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "C"] <- "S"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "G"] <- "C"
snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "S"] <- "G"

snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "A"] <- "S"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "T"] <- "A"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "S"] <- "T"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "C"] <- "S"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "G"] <- "C"
snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "S"] <- "G"


snv_info$GSM2 <- rep("no-change", length(snv_info$GenomonSplicingMutation))
snv_info$GSM2[snv_info$GenomonSplicingMutation != "---"] <- snv_info$GenomonSplicingMutation[snv_info$GenomonSplicingMutation != "---"]
snv_info$GSM2[grep(";", snv_info$GenomonSplicingMutation)] <- "complex"
snv_info$GSM2[snv_info$GSM2 == "opposite-side-intron-retention"] <- "intron-retention"
snv_info$GSM2[snv_info$GSM2 == "intronic-alternative-5'-splice-site"] <- "alternative-5'-splice-site"
snv_info$GSM2[snv_info$GSM2 == "intronic-alternative-3'-splice-site"] <- "alternative-3'-splice-site"

# snv_info$GSM2[snv_info$GSM2 == "alternative-5'-splice-site"] <- "alternative-splice-site"
# snv_info$GSM2[snv_info$GSM2 == "alternative-3'-splice-site"] <- "alternative-splice-site"


##########
# donor maxEnt score
snv_info_d <- snv_info %>% filter(Type_Motif == "donor")

ent_chr_d <- paste("chr", snv_info_d$Chr_Motif, sep = "")

ent_start_d <- rep(0, length(snv_info_d$Start_Motif))
ent_start_d[snv_info_d$Strand_Motif == "+"] <- snv_info_d$Start_Motif[snv_info_d$Strand_Motif == "+"] + exon_size_d - 3
ent_start_d[snv_info_d$Strand_Motif == "-"] <- snv_info_d$Start_Motif[snv_info_d$Strand_Motif == "-"] + intron_size_d - 8

ent_end_d <- rep(0, length(snv_info_d$Start_Motif))
ent_end_d[snv_info_d$Strand_Motif == "+"] <- snv_info_d$End_Motif[snv_info_d$Strand_Motif == "+"] - intron_size_d + 8
ent_end_d[snv_info_d$Strand_Motif == "-"] <- snv_info_d$End_Motif[snv_info_d$Strand_Motif == "-"]- exon_size_d + 3


gr <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = ent_chr_d, start = ent_start_d, end = ent_end_d, 
             strand = snv_info_d$Strand_Motif), ignore.strand = FALSE)

donor_seq_wt <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr)

donor_mes_wt <- unlist(lapply(as.character(donor_seq_wt), function(x) {score5(mes, x, pos = 3)}))

hb_score_wt <- hbond(hb, as.character(donor_seq_wt), 3)


mut_pos_d <- snv_info_d$Rel_Start_Motif + 3 - exon_size_d
mut_alt_d <- snv_info_d$Alt_Mut

donor_seq_mt <- unlist(lapply(1:length(donor_seq_wt), 
                              function(x) {
                                paste(substring(donor_seq_wt[x], 1, mut_pos_d[x] - 1), 
                                      mut_alt_d[x], 
                                      substring(donor_seq_wt[x], mut_pos_d[x] + 1, 11), 
                                      sep ="")
                              }))

donor_mes_mt <- unlist(lapply(as.character(donor_seq_mt), function(x) {score5(mes, x, pos = 3)}))
hb_score_mt <- hbond(hb, as.character(donor_seq_mt), 3)
##########

##########
# acceptor masEntScore
snv_info_a <- snv_info %>% filter(Type_Motif == "acceptor")

ent_chr_a <- paste("chr", snv_info_a$Chr_Motif, sep = "")

ent_start_a <- rep(0, length(snv_info_a$Start_Motif))
ent_start_a[snv_info_a$Strand_Motif == "+"] <- snv_info_a$Start_Motif[snv_info_a$Strand_Motif == "+"] + intron_size_a - 20 
ent_start_a[snv_info_a$Strand_Motif == "-"] <- snv_info_a$Start_Motif[snv_info_a$Strand_Motif == "-"] + exon_size_a - 3

ent_end_a <- rep(0, length(snv_info_a$Start_Motif))
ent_end_a[snv_info_a$Strand_Motif == "+"] <- snv_info_a$End_Motif[snv_info_a$Strand_Motif == "+"] - exon_size_a + 3
ent_end_a[snv_info_a$Strand_Motif == "-"] <- snv_info_a$End_Motif[snv_info_a$Strand_Motif == "-"] - intron_size_a + 20 


gr <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = ent_chr_a, start = ent_start_a, end = ent_end_a, 
             strand = snv_info_a$Strand_Motif), ignore.strand = FALSE)


acceptor_seq_wt <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr)

acceptor_mes_wt <- unlist(lapply(as.character(acceptor_seq_wt), function(x) {score3(mes, x, pos = 20)}))


mut_pos_a <- snv_info_a$Rel_Start_Motif + 20 - intron_size_a
mut_alt_a <- snv_info_a$Alt_Mut

acceptor_seq_mt <- unlist(lapply(1:length(acceptor_seq_wt), 
                              function(x) {
                                paste(substring(acceptor_seq_wt[x], 1, mut_pos_a[x] - 1), 
                                      mut_alt_a[x], 
                                      substring(acceptor_seq_wt[x], mut_pos_a[x] + 1, 23), 
                                      sep ="")
                              }))

acceptor_mes_mt <- unlist(lapply(as.character(acceptor_seq_mt), function(x) {score3(mes, x, pos = 20)}))
##########


is_gsm <- factor(c(ifelse(snv_info_d$GSM2 != "no-change", "Aberrant", "Normal"), 
                   ifelse(snv_info_a$GSM2 != "no-change", "Aberrant", "Normal")),
                 levels = c("Aberrant", "Normal"))



mes_df <- data.frame(motif_type = factor(c(rep("donor", nrow(snv_info_d)), rep("acceptor", nrow(snv_info_a))),
                                         levels = c("donor", "acceptor"),
                                         labels = c("Donor", "Acceptor")),
                     mes_diff = c(donor_mes_mt - donor_mes_wt, acceptor_mes_mt - acceptor_mes_wt),
                     splice_class = factor(c(snv_info_d$GSM2, snv_info_a$GSM2),
                                           levels = rev(c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
                                                      "intron-retention", "complex", "no-change")),
                                           labels = rev(c("Exon skip", "Alternative 5'-ss", "Alternative 3'-ss",
                                                          "Intron retention", "Complex", "No change"))),
                     mut_pos = c(snv_info_d$Rel_Start_Motif, snv_info_a$Rel_Start_Motif),
                     is_gsm = is_gsm)

g_mes_d <- ggplot(mes_df %>% filter(splice_class != "Alternative 3'-ss" & motif_type == "Donor"), aes(x = splice_class, y = mes_diff, fill = splice_class)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  coord_flip() +
  ggtitle("Donor disruption") +
  ylim(c(-15, 10)) +
  my_theme() +
  scale_fill_manual(values = splicing_class_colour) + 
  labs(x = "", y = "") +
  guides(fill = FALSE)

g_mes_a <- ggplot(mes_df %>% filter(splice_class != "Alternative 5'-ss" & motif_type == "Acceptor"), aes(x = splice_class, y = mes_diff, fill = splice_class)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  coord_flip() +
  ggtitle("Acceptor disruption") +
  labs(x = "", y = "") +
  ylim(c(-15, 10)) +
  my_theme() +
  scale_fill_manual(values = splicing_class_colour) + 
  labs(x = "", y = "") +
  guides(fill = FALSE)

ylabel <- ggdraw() + draw_label("Diff. of MaxEnt score")

plot_grid(g_mes_d, g_mes_a, ylabel, ncol = 1, rel_heights = c(1, 1, 0.1))

ggsave("../figure/diff_mes_spliceclass.pdf", width = 4, height = 4)


pos_colour <- rep("grey30", 9)
pos_colour[4:5] <- "red"

ggplot(mes_df %>% filter(motif_type == "Donor"), 
       aes(x = factor(mut_pos, levels = 1:9), y = mes_diff, fill = is_gsm)) +
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  labs(fill = "") +
  my_theme() +
  ylim(c(-15, 10)) + 
  labs(x = "Intronic position", y = "Diff. of MaxEnt score", fill = "Splicing") +
  scale_x_discrete(limits = 1:9, 
                   labels = c("-3", "-2", "-1", "1", "2", "3", "4", "5", "6")) +
  scale_fill_manual(values = c("#ef8a62", "#999999")) +
  theme(axis.text.x = element_text(colour = pos_colour),
        legend.position = "bottom")

ggsave("../figure/diff_mes_mutpos_donor.pdf", width = 6, height = 4)



pos_colour <- rep("grey30", 7)
pos_colour[5:6] <- "red"

ggplot(mes_df %>% filter(motif_type == "Acceptor"), 
       aes(x = factor(mut_pos, levels = 1:7), y = mes_diff, fill = is_gsm)) +
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  my_theme() +
  labs(fill = "") +
  ylim(c(-15, 10)) + 
  labs(x = "Intronic position", y = "Diff. of MaxEnt score", fill = "Splicing") +
  scale_x_discrete(limits = 1:7, 
                   labels = c("6", "5", "4", "3", "2", "1", "-1")) +
  scale_fill_manual(values = c("#ef8a62", "#999999")) +
  theme(axis.text.x = element_text(colour = pos_colour),
        legend.position = "bottom") 

ggsave("../figure/diff_mes_mutpos_acceptor.pdf", width = 6, height = 4)



is_gsm <- factor(ifelse(snv_info_d$GSM2 != "no-change", "Aberrant", "Normal"),
                 levels = c("Aberrant", "Normal"))

hb_df <- data.frame(hb_diff = c(hb_score_mt - hb_score_wt),
                    splice_class = factor(snv_info_d$GSM2,
                                          levels = rev(c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
                                                         "intron-retention", "complex", "no-change")),
                                          labels = rev(c("Exon skip", "Alternative 5'-ss", "Alternative 3'-ss",
                                                         "Intron retention", "Complex", "No change"))),
                     mut_pos = snv_info_d$Rel_Start_Motif,
                     is_gsm = is_gsm)

ggplot(hb_df %>% filter(splice_class != "Alternative 3' splice site"), aes(x = splice_class, y = hb_diff, fill = splice_class)) + 
  geom_boxplot(size = 0.3, outlier.size = 0.4) +
  coord_flip() +
  ylim(c(-20, 10)) +
  my_theme() +
  scale_fill_manual(values = splicing_class_colour) + 
  labs(x = "", y = "Diff. of HBond score") +
  guides(fill = FALSE)


ggsave("../figure/diff_hb_spliceclass.pdf", width = 6, height = 2.4)





pos_colour <- rep("grey30", 10)
pos_colour[4:5] <- "red"

ggplot(hb_df, 
       aes(x = factor(mut_pos, levels = 1:9), y = hb_diff, fill = is_gsm)) +
  geom_boxplot(outlier.size = 0.4) +
  my_theme() +
  labs(fill = "") +
  ylim(c(-20, 10)) + 
  labs(x = "Intronic position", y = "Diff. of HBond score", fill = "Splicing") +
  scale_x_discrete(limits = 1:9, 
                   labels = c("-3", "-2", "-1", "1", "2", "3", "4", "5", "6")) +
  theme(axis.text.x = element_text(colour = pos_colour),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#ef8a62", "#999999")) 

ggsave("../figure/diff_hb_mutpos_donor.pdf", width = 6, height = 4)



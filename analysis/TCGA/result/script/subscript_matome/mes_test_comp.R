library(dplyr)
library(spliceSites)
library(ggplot2)
library(cowplot)

source("../../../conf/plot_config.R")

mes<-load.maxEnt()
hb<-load.hbond()


get_print_info <- function(score_df, is_mes, is_donor, is_dummy) {

  if (is_donor == TRUE) {
    ind2sptype <- c("Exon skipping", "Alternative 5'SS",
                    "Intron retention", "Complex", "Normal splicing")
  } else {
    ind2sptype <- c("Exon skipping", "Alternative 3'SS",
                    "Intron retention", "Complex", "Normal splicing")
  }

  tmotif <- ifelse(is_donor, "donor", "acceptor")

  Ps <- c()


  for(i1 in 1:4) {
    for (i2 in (i1 + 1):5) {

      score_1 <-  score_df %>% filter(Type_Motif == tmotif, Splice_Class == ind2sptype[i1])
      score_2 <-  score_df %>% filter(Type_Motif == tmotif, Splice_Class == ind2sptype[i2])

      if (is_mes == TRUE) {

        W_g <- wilcox.test(score_1$mes, score_2$mes, alternative = "greater")
        W_l <- wilcox.test(score_1$mes, score_2$mes, alternative = "less")
        # mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], mlogp = mlogp, motif = tmotif))

      } else {

        W_g <- wilcox.test(score_1$hbs, score_2$hbs, alternative = "greater")
        W_l <- wilcox.test(score_1$hbs, score_2$hbs, alternative = "less")
        # mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], mlogp = mlogp, motif = tmotif))

      }
    }
  }

  Ps$mlogp <- pmin(Ps$mlogp, 10)
  Ps$mlogp <- pmax(Ps$mlogp, -10)
  # Ps[Ps$x == Ps$y,"mlogp"] <- 0

  print(Ps)

  if (is_mes == TRUE) {
    if (is_donor == TRUE) {
      ttitle <- "MaxExt score, donor"
    } else {
      ttitle <- "MaxEnt score, acceptor"
    }
  } else {
    ttitle <- "H-bond score, donor"
  }


  p <- ggplot(Ps, aes(x = x, y = y, fill = mlogp)) +
    geom_tile(colour = "grey30") +
    ggtitle(ttitle) +
    coord_flip() +
    # facet_grid(.~type) +
    my_theme() +
    scale_fill_gradient2(low = "#2166ac", mid = "#ffffff", high = "#b2182b") +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.key.width = unit(1, "cm"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          plot.background = element_blank(),
          legend.position = "bottom")


  # labs(x = "Splicing pattern 1", y = "Splicing pattern 2", fill = "Differential score for splicing pattern 1 v.s. 2")

  if (is_dummy == FALSE) {
    p <- p + guides(fill = FALSE) +
      labs(x = "Splicing pattern 1", y = "Splicing pattern 2")
  } else {
    p <- p + labs(fill = "Differential score for splicing pattern 1 v.s. 2")
  }

  return(p)
}




snv_info <- read.table("../temporary/TCGA.savnet.allele_count.summary.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "") %>%
  filter(Ref_Mut != "-" & Alt_Mut != "-")

# snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "A"] <- "S"
# snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "T"] <- "A"
# snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "S"] <- "T"
# snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "C"] <- "S"
# snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "G"] <- "C"
# snv_info$Ref_Mut[snv_info$Strand_Motif == "-" & snv_info$Ref_Mut == "S"] <- "G"

# snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "A"] <- "S"
# snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "T"] <- "A"
# snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "S"] <- "T"
# snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "C"] <- "S"
# snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "G"] <- "C"
# snv_info$Alt_Mut[snv_info$Strand_Motif == "-" & snv_info$Alt_Mut == "S"] <- "G"



snv_info$GSM2 <- rep("no-change", length(snv_info$GenomonSplicingMutation))
snv_info$GSM2[snv_info$GenomonSplicingMutation != "---"] <- snv_info$GenomonSplicingMutation[snv_info$GenomonSplicingMutation != "---"]
snv_info$GSM2[grep(";", snv_info$GenomonSplicingMutation)] <- "complex"
snv_info$GSM2[snv_info$GSM2 == "opposite-side-intron-retention"] <- "intron-retention"
snv_info$GSM2[snv_info$GSM2 == "intronic-alternative-5'-splice-site"] <- "alternative-5'-splice-site"
snv_info$GSM2[snv_info$GSM2 == "intronic-alternative-3'-splice-site"] <- "alternative-3'-splice-site"

# snv_info$GSM2[snv_info$GSM2 == "alternative-5'-splice-site"] <- "alternative-splice-site"
# snv_info$GSM2[snv_info$GSM2 == "alternative-3'-splice-site"] <- "alternative-splice-site"

snv_info <- snv_info %>% filter(GSM2 != "no-change" | FPKM >= 10)

print(nrow(snv_info))

snv_info$Splice_Class <- factor(snv_info$GSM2,
                          levels = c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
                                     "intron-retention", "complex", "no-change"),
                          labels = c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                                     "Intron retention", "Complex", "Normal splicing"))


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

snv_info_d$mes <- donor_mes_wt

hb_score_wt <- hbond(hb, as.character(donor_seq_wt), 3)

snv_info_d$hbs <- hb_score_wt


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

snv_info_a$mes <- acceptor_mes_wt

##########


##########
# for legend

tdf <- data.frame(x = c(0, 1), y = c(0, 1), z = c(0, 10))

p_legend_less <- ggplot(tdf, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(high = "#2166ac", low = "#ffffff") +
  my_theme() +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.6, "cm")) +
  labs(fill = "-Log10(P-value) (Splicing patttern 1 < 2)")

p_legend_greater <- ggplot(tdf, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(high = "#b2182b", low = "#ffffff") +
  my_theme() +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.6, "cm")) +
  labs(fill = "-Log10(P-value) (Splicing patttern 1 > 2)")
##########

p_mes_donor <- get_print_info(snv_info_d, TRUE, TRUE, FALSE)
p_hb_donor <- get_print_info(snv_info_d, FALSE, TRUE, FALSE)
p_mes_acceptor <- get_print_info(snv_info_a, TRUE, FALSE, FALSE)

plot_grid(plot_grid(p_mes_donor, p_mes_acceptor, p_hb_donor, nrow = 1),
          plot_grid(g_legend(p_legend_less), g_legend(p_legend_greater), nrow = 1),
          ncol = 1, rel_heights = c(1.0, 0.2))

ggsave("../figure/mes_hb_wt_comp.tiff, width = 18, height = 7.5, dpi = 600, units = "cm")



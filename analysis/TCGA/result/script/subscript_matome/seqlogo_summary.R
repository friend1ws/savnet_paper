library(dplyr)
# library(seqLogo)
library(ggplot2)
library(ggseqlogo)
library(cowplot)
library(Cairo)

Cairo()

source("../../../conf/plot_config.R")

donor_pos_label <- c(add_emdash("3"), add_emdash("2"), add_emdash("1"), "+1", "+2", "+3", "+4", "+5", "+6")
acceptor_pos_label <- c("+6", "+5", "+4", "+3", "+2", "+1", add_emdash("1"))


splicing_mutation <- read.table("../temporary/TCGA.savnet.allele_count.summary.txt", sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE) %>%
  filter(Ref_Mut != "-" & Alt_Mut != "-") 


splicing_mutation$GSM2 <- splicing_mutation$GenomonSplicingMutation
splicing_mutation$GSM2[splicing_mutation$GSM2 == "---"] <- "no-change"
splicing_mutation$GSM2[splicing_mutation$GSM2 %in% c("intron-retention", "opposite-side-intron-retention",
                                                     "intron-retention;opposite-side-intron-retention",
                                                     "opposite-side-intron-retention;intron-retention")] <- "intron-retention"
splicing_mutation$GSM2[grep(";", splicing_mutation$GSM2)] <- "complex"
splicing_mutation$GSM2[splicing_mutation$GSM2 == "intronic-alternative-5'-splice-site"] <- "alternative-5'-splice-site"
splicing_mutation$GSM2[splicing_mutation$GSM2 == "intronic-alternative-3'-splice-site"] <- "alternative-3'-splice-site"

splicing_mutation$GSM2 <- factor(splicing_mutation$GSM2,
                                 levels = rev(c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
                                                "intron-retention", "complex", "no-change")),
                                 labels = rev(c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                                                "Intron retention", "Complex", "Normal splicing"))
)


theme_nonbottom <- function() {
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "lines"),
          axis.text = element_text(size = 7),
          # axis.title = element_text(size = 7),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title = element_blank())
}


theme_bottom <- function() {
    theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "lines"),
          panel.border = element_blank(),
          axis.text = element_text(size = 7),
          # axis.title = element_text(size = 7),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4)) # +
     # scale_x_continuous(breaks = 1:9,
     #                    labels = c("-3", "-2", "-1", "1", "2", "3", "4", "5", "6"))

}

# print(theme_nonbottom())
# print("OK")
# print(theme_bottom())
# print("OK2")

get_base_count_mat <- function(motif_seqs, nbase) {

  base_count_mat <- matrix(0, nbase, 4)
  for(n in 1:nbase) {
    base_count <- table(substring(motif_sub$Motif_Seq, n + 0, n + 0))
    base_count_mat[n, 1] <- base_count["A"]
    base_count_mat[n, 2] <- base_count["C"]
    base_count_mat[n, 3] <- base_count["G"]
    base_count_mat[n, 4] <- base_count["T"]
  }
  base_count_mat[is.na(base_count_mat)] <- 0  
  return(base_count_mat)
}




seq_logo_print_list <- list(plot_grid(ggdraw() + draw_label(""),
                                      ggdraw() + draw_label(""),
                                      ggdraw() + draw_label("Normal splicing", size = 7),
                                      ggdraw() + draw_label("Abnormal splicing", size = 7),
                                      ncol = 4, rel_widths = c(0.15, 0.15, 1, 1)))


for(pos in c(2, 3, 7, 8, 9)) {
  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Normal splicing") %>% filter(FPKM >= 10.0)
  
  tgg_logo_no <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
    labs(x = "", y = "Bits") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))

  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 != "Normal splicing")
  
  tgg_logo_yes <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
    labs(x = "", y = "Bits") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))


  if (pos != 9) {
        tgg_logo_no <- tgg_logo_no + theme_nonbottom()
        tgg_logo_yes <- tgg_logo_yes + theme_nonbottom()
  } else {
        tgg_logo_no <- tgg_logo_no + theme_bottom() + 
            scale_x_continuous(breaks = 1:9, labels = donor_pos_label)
        tgg_logo_yes <- tgg_logo_yes + theme_bottom() + 
            scale_x_continuous(breaks = 1:9, labels = donor_pos_label)
  }

  label_pos <- ifelse(pos > 3, paste("+", as.character(pos - 3), sep = ""), add_emdash(as.character(4 - pos)))
  bits_label <- ifelse(pos == 7, "Bits", "")
  seq_logo_print_list <- c(seq_logo_print_list, list(
    plot_grid(ggdraw() + draw_label(label_pos, size = 7), 
              ggdraw() + draw_label(bits_label, angle = 90, size = 7),
              tgg_logo_no, 
              tgg_logo_yes, ncol = 4, rel_widths = c(0.15, 0.15, 1, 1))
  )
  )
  
}


# seq_logo_print_list <- c(seq_logo_print_list, list(
#   plot_grid(ggdraw() + draw_label("", size = 7),
#             ggdraw() + draw_label("", size = 7),
#             ggdraw() + draw_label("Splicing donor motif", size = 7, hjust = 0.5),
#             rel_widths = c(0.15, 0.15, 2))
# )
# )           

# plot_grid(ggdraw() + draw_label("Variant position at splicing donor sites", angle = 90, size = 7),
#           plot_grid(plotlist = seq_logo_print_list, ncol = 1, rel_heights = c(0.4, 1, 1, 1, 1, 1.2, 0.25)),
#           ncol = 2, rel_width = c(0.06, 1))

plot_grid(ggdraw() + draw_label("Variant position at splicing donor sites", angle = 90, size = 7),
          plot_grid(
            plot_grid(plotlist = seq_logo_print_list, ncol = 1, rel_heights = c(0.4, 1, 1, 1, 1, 1.2)),
            ggdraw() + draw_label("                      Splicing donor motif", size = 7), 
            align = "v", ncol = 1, rel_heights = c(5.7, 0.25)),
            ncol = 2, rel_widths = c(0.06, 1))


ggsave("../figure/seqlogo_list_simple.tiff", width = 10, height = 9.5, dpi = 600, units = "cm")


##########
# detailed version

count_thres <- 25 
seq_logo_print_list <- list(plot_grid(ggdraw() + draw_label(""),
                                      ggdraw() + draw_label(""),
                                      ggdraw() + draw_label("Normal splicing", size = 7),
                                      ggdraw() + draw_label("Exon skipping", size = 7),
                                      ggdraw() + draw_label("Alternative 5'SS", size = 7),
                                      ggdraw() + draw_label("Intron retention", size = 7),
                                      ggdraw() + draw_label("Complex", size = 7),
                                      ncol = 7, rel_widths = c(0.15, 0.15, 1, 1, 1, 1, 1)))

count_mat <- matrix(0, 9, 5)

for(pos in c(2, 3, 7, 8, 9)) {
  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Normal splicing") %>% filter(FPKM >= 10.0)
  count_mat[pos, 1] <- nrow(motif_sub)

  if (count_mat[pos, 1] >= count_thres) {
    tgg_logo_no <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") 
  } else {
    tgg_logo_no <- ggplot() 
  }

  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Exon skipping")
  count_mat[pos, 2] <- nrow(motif_sub)
  
  if (count_mat[pos, 2] >= count_thres) {
    tgg_logo_ES <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") 
  } else {
    tgg_logo_ES <- ggplot()
  }

  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Alternative 5'SS")
  count_mat[pos, 3] <- nrow(motif_sub)
  
  if (count_mat[pos, 3] >= count_thres) {
    tgg_logo_A5S <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "")
  } else {
    tgg_logo_A5S <- ggplot()
  }
  

  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Intron retention")
  count_mat[pos, 4] <- nrow(motif_sub)
 
  if (count_mat[pos, 4] >= count_thres) {
    tgg_logo_IR <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "")
  } else {
    tgg_logo_IR <- ggplot()
  }
  

  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Complex")
  count_mat[pos, 5] <- nrow(motif_sub)


  if (count_mat[pos, 5] >= count_thres) {
    tgg_logo_CP <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") 
  } else {
    tgg_logo_CP <- ggplot() 
  }
  
  if (pos != 9) {
      tgg_logo_no <- tgg_logo_no + theme_nonbottom() + scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      tgg_logo_ES <- tgg_logo_ES + theme_nonbottom() + scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      tgg_logo_A5S <- tgg_logo_A5S + theme_nonbottom() + scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      tgg_logo_IR <- tgg_logo_IR + theme_nonbottom() + scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      tgg_logo_CP <- tgg_logo_CP + theme_nonbottom() + scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
  } else {
      if (count_mat[pos, 1] >= count_thres) tgg_logo_no <- tgg_logo_no + theme_bottom() + 
                                              scale_x_continuous(breaks = 1:9, labels = donor_pos_label) + 
                                              scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      if (count_mat[pos, 2] >= count_thres) tgg_logo_ES <- tgg_logo_ES + theme_bottom() + 
                                              scale_x_continuous(breaks = 1:9, labels = donor_pos_label) + 
                                              scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      if (count_mat[pos, 3] >= count_thres) tgg_logo_A5S <- tgg_logo_A5S + theme_bottom() + 
                                              scale_x_continuous(breaks = 1:9, labels = donor_pos_label) +
                                              scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      if (count_mat[pos, 4] >= count_thres) tgg_logo_IR <- tgg_logo_IR + theme_bottom() + 
                                              scale_x_continuous(breaks = 1:9, labels = donor_pos_label) +
                                              scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))
      if (count_mat[pos, 5] >= count_thres) tgg_logo_CP <- tgg_logo_CP + theme_bottom() + 
                                              scale_x_continuous(breaks = 1:9, labels = donor_pos_label) +
                                              scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))  
  }


 
  label_pos <- ifelse(pos > 3, paste("+", as.character(pos - 3), sep = ""), add_emdash(as.character(4 - pos)))
  bits_label <- ifelse(pos == 7, "Bits", "")
  seq_logo_print_list <- c(seq_logo_print_list, list(
    plot_grid(ggdraw() + draw_label(label_pos, size = 7), 
              ggdraw() + draw_label(bits_label, angle = 90, size = 7),
              tgg_logo_no, 
              tgg_logo_ES,
              tgg_logo_A5S,
              tgg_logo_IR,
              tgg_logo_CP, 
              ncol = 7, rel_widths = c(0.15, 0.15, 1, 1, 1, 1, 1))
    )
  )
  
}

# axis_x <- ggplot(data.frame(x = c(1, 9), y = c(0, 0)), 
#                  aes(x = x, y = y)) + geom_polygon() +
#   labs(x = "") +
#   scale_x_continuous(breaks = 1:9,
#                      labels = c("-3", "-2", "-1", "1", "2", "3", "4", "5", "6")) +
#   theme(axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank())


# seq_logo_print_list <- c(seq_logo_print_list, list(
#   plot_grid(ggdraw() + draw_label(""), 
#             axis_x, axis_x, axis_x, axis_x, axis_x, ncol = 6, rel_widths = c(0.15, 1, 1, 1, 1, 1))
# )
# )


plot_grid(ggdraw() + draw_label("Variant position at splicing donor sites", angle = 90, size = 7),
          plot_grid(
            plot_grid(plotlist = seq_logo_print_list, ncol = 1, rel_heights = c(0.4, 1, 1, 1, 1, 1.2)),
            ggdraw() + draw_label("                    Splicing donor motif", size = 7), 
            ncol = 1, rel_heights = c(5.7, 0.25)),
          ncol = 2, rel_widths = c(0.025, 1))


ggsave("../figure/seqlogo_list_detail.tiff", width = 20, height = 9.5, dpi = 600, units = "cm")




seq_logo_print_list <- list(plot_grid(ggdraw() + draw_label(""),
                                      ggdraw() + draw_label(""),
                                      ggdraw() + draw_label("Normal splicing", size = 7),
                                      ggdraw() + draw_label("Abnormal splicing", size = 7),
                                      ncol = 4, rel_widths = c(0.15, 0.15, 1, 1)))

for(pos in c(1, 2, 3, 4, 7)) {

  motif_sub <- splicing_mutation %>%
    filter(Type_Motif == "acceptor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Normal splicing") %>% filter(FPKM >= 10.0)

  tgg_logo_no <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 7)), TRUE) +
    labs(x = "", y = "") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))

  motif_sub <- splicing_mutation %>%
    filter(Type_Motif == "acceptor", Rel_Start_Motif == pos) %>%
    filter(GSM2 != "Normal splicing")

  tgg_logo_yes <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 7)), TRUE) +
    labs(x = "", y = "") +
    scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0))


  if (pos != 7) {
        tgg_logo_no <- tgg_logo_no + theme_nonbottom()
        tgg_logo_yes <- tgg_logo_yes + theme_nonbottom()
  } else {
        tgg_logo_no <- tgg_logo_no + theme_bottom() + scale_x_continuous(breaks = 1:7, labels = acceptor_pos_label)
        tgg_logo_yes <- tgg_logo_yes + theme_bottom() + scale_x_continuous(breaks = 1:7, labels = acceptor_pos_label)
  }

  label_pos <- ifelse(pos > 6, add_emdash(as.character(pos - 6)), paste("+", as.character(-pos + 7), sep = ""))
  bits_label <- ifelse(pos == 3, "Bits", "")
  seq_logo_print_list <- c(seq_logo_print_list, list(
    plot_grid(ggdraw() + draw_label(label_pos, size = 7),
              ggdraw() + draw_label(bits_label, angle = 90, size = 7),
              tgg_logo_no,
              tgg_logo_yes, ncol = 4, rel_widths = c(0.15, 0.15, 0.9, 0.9))
  )
  )

}


plot_grid(ggdraw() + draw_label("Variant position at splicing acceptor sites", angle = 90, size = 7),
          plot_grid(
            plot_grid(plotlist = seq_logo_print_list, ncol = 1, rel_heights = c(0.4, 1, 1, 1, 1, 1.2)),
            ggdraw() + draw_label("                    Splicing acceptor motif", size = 7),
            align = "v", ncol = 1, rel_heights = c(5.7, 0.25)),
          ncol = 2, rel_widths = c(0.06, 1))


ggsave("../figure/seqlogo_acceptor.tiff", width = 10, height = 9.5, dpi = 600, units = "cm")


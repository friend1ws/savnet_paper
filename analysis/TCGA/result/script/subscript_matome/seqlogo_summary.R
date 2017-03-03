library(dplyr)
library(seqLogo)
library(ggplot2)
library(ggseqlogo)
library(cowplot)

# source("ggseqlogo.R")

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
                                 labels = rev(c("Exon skip", "Alternative 5'-ss", "Alternative 3'-ss",
                                                "Intron retention", "Complex", "No change"))
)

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
                                      ggdraw() + draw_label("No change"),
                                      ggdraw() + draw_label("Abnormal splicing"),
                                      ncol = 3, rel_widths = c(0.15, 1, 1)))


for(pos in c(2, 3, 7, 8, 9)) {
  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "No change") %>% filter(FPKM >= 10.0)
  
  tgg_logo_no <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
    labs(x = "", y = "") +
    theme_nothing() 

                           
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 != "No change")
  
  tgg_logo_yes <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
    labs(x = "", y = "") +
    theme_nothing() 

  label_pos <- ifelse(pos > 3, as.character(pos - 3), as.character(pos - 4))
  seq_logo_print_list <- c(seq_logo_print_list, list(
    plot_grid(ggdraw() + draw_label(label_pos), 
              tgg_logo_no, 
              tgg_logo_yes, ncol = 3, rel_widths = c(0.15, 1, 1))
  )
  )
  
}

axis_x <- ggplot(data.frame(x = c(1, 9), y = c(0, 0)), 
                 aes(x = x, y = y)) + geom_polygon() +
  labs(x = "") +
  scale_x_continuous(breaks = 1:9,
                     labels = c("-3", "-2", "-1", "1", "2", "3", "4", "5", "6")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


seq_logo_print_list <- c(seq_logo_print_list, list(
  plot_grid(ggdraw() + draw_label(""), axis_x, axis_x, ncol = 3, rel_widths = c(0.15, 1, 1))
)
)


plot_grid(ggdraw() + draw_label("Mutation position at splicing donor sites", angle = 90),
          plot_grid(
            plot_grid(plotlist = seq_logo_print_list, ncol = 1, rel_heights = c(0.4, 1, 1, 1, 1, 1, 0.3)),
            ggdraw() + draw_label("Position"), 
            ncol = 1, rel_heights = c(5.7, 0.25)),
          ncol = 2, rel_widths = c(0.05, 1))


ggsave("../figure/seqlogo_list_simple.pdf", width = 6, height = 6)




##########
# detailed version

count_thres <- 30
seq_logo_print_list <- list(plot_grid(ggdraw() + draw_label(""),
                                      ggdraw() + draw_label("No change"),
                                      ggdraw() + draw_label("Exon skip"),
                                      ggdraw() + draw_label("Alternative 5'-ss"),
                                      ggdraw() + draw_label("Intron retention"),
                                      ggdraw() + draw_label("Complex"),
                                      ncol = 6, rel_widths = c(0.15, 1, 1, 1, 1, 1)))

count_mat <- matrix(0, 9, 5)

for(pos in c(2, 3, 7, 8, 9)) {
  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "No change") %>% filter(FPKM >= 10.0)
  count_mat[pos, 1] <- nrow(motif_sub)

  if (count_mat[pos, 1] >= count_thres) {
    tgg_logo_no <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") +
      theme_nothing()  
  } else {
    tgg_logo_no <- ggplot() + theme_nothing()
  }
  
  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Exon skip")
  count_mat[pos, 2] <- nrow(motif_sub)
  
  if (count_mat[pos, 2] >= count_thres) {
    tgg_logo_ES <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") +
      theme_nothing()  
  } else {
    tgg_logo_ES <- ggplot() + theme_nothing()
  }

  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Alternative 5'-ss")
  count_mat[pos, 3] <- nrow(motif_sub)
  
  if (count_mat[pos, 3] >= count_thres) {
    tgg_logo_A5S <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") +
      theme_nothing()  
  } else {
    tgg_logo_A5S <- ggplot() + theme_nothing()
  }
  

  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Intron retention")
  count_mat[pos, 4] <- nrow(motif_sub)
 
  if (count_mat[pos, 4] >= count_thres) {
    tgg_logo_IR <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") +
      theme_nothing()  
  } else {
    tgg_logo_IR <- ggplot() + theme_nothing()
  }
  
 
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GSM2 == "Complex")
  count_mat[pos, 5] <- nrow(motif_sub)

  if (count_mat[pos, 5] >= count_thres) {
    tgg_logo_CP <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
      labs(x = "", y = "") +
      theme_nothing()  
  } else {
    tgg_logo_CP <- ggplot() + theme_nothing()
  }
  
  tgg_logo_CP <- ggseqlogo(t(get_base_count_mat(motif_sub$Motif_Seq, 9)), TRUE) + 
    labs(x = "", y = "") +
    theme_nothing()   
  
 
  label_pos <- ifelse(pos > 3, as.character(pos - 3), as.character(pos - 4))
  seq_logo_print_list <- c(seq_logo_print_list, list(
    plot_grid(ggdraw() + draw_label(label_pos), 
              tgg_logo_no, 
              tgg_logo_ES,
              tgg_logo_A5S,
              tgg_logo_IR,
              tgg_logo_CP, 
              ncol = 6, rel_widths = c(0.15, 1, 1, 1, 1, 1))
  )
  )
  
}

axis_x <- ggplot(data.frame(x = c(1, 9), y = c(0, 0)), 
                 aes(x = x, y = y)) + geom_polygon() +
  labs(x = "") +
  scale_x_continuous(breaks = 1:9,
                     labels = c("-3", "-2", "-1", "1", "2", "3", "4", "5", "6")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


seq_logo_print_list <- c(seq_logo_print_list, list(
  plot_grid(ggdraw() + draw_label(""), 
            axis_x, axis_x, axis_x, axis_x, axis_x, ncol = 6, rel_widths = c(0.15, 1, 1, 1, 1, 1))
)
)


plot_grid(ggdraw() + draw_label("Mutation position at splicing donor sites", angle = 90),
          plot_grid(
            plot_grid(plotlist = seq_logo_print_list, ncol = 1, rel_heights = c(0.4, 1, 1, 1, 1, 1, 0.3)),
            ggdraw() + draw_label("Position"), 
            ncol = 1, rel_heights = c(5.7, 0.25)),
          ncol = 2, rel_widths = c(0.03, 1))


ggsave("../figure/seqlogo_list_detail.pdf", width = 11, height = 6)




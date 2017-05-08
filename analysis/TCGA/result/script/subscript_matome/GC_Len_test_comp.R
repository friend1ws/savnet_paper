library(ggplot2)
library(dplyr)
library(cowplot)

source("../../../conf/plot_config.R")

GC <- read.table("../temporary/TCGA.savnet.gc_intron.txt", sep = "\t", header = TRUE, quote = "")

GC <- GC %>% filter(Splice_Class != "no-change" | FPKM >= 10.0)

print(nrow(GC))

GC$Splice_Class[GC$Splice_Class == "intronic-alternative-5'-splice-site"] <- "alternative-5'-splice-site"
GC$Splice_Class[GC$Splice_Class == "intronic-alternative-3'-splice-site"] <- "alternative-3'-splice-site"
GC$Splice_Class[GC$Splice_Class == "opposite-side-intron-retention"] <- "intron-retention"


GC$Splice_Class <- factor(GC$Splice_Class,
                          levels = c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site", 
                                     "intron-retention", "complex", "no-change"),
                          labels = c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
                                     "Intron retention", "Complex", "Normal splicing"))


ind2sptype <- c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS",
               "Intron retention", "Complex", "Normal splicing")
# ind2motif <- c("donor", "acceptor")


get_print_info <- function(GC_Len_df, is_donor, is_GC, is_dummy) {

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
      
      GC_1 <-  GC_Len_df %>% filter(Type_Motif == tmotif, Splice_Class == ind2sptype[i1])
      GC_2 <-  GC_Len_df %>% filter(Type_Motif == tmotif, Splice_Class == ind2sptype[i2])
      
      if (is_GC == TRUE) {
        
        W_g <- wilcox.test(GC_1$GC_intron_5prime, GC_2$GC_intron_5prime, alternative = "greater")
        W_l <- wilcox.test(GC_1$GC_intron_5prime, GC_2$GC_intron_5prime, alternative = "less")
        mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        # mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], type = "GC_intron_5prime", mlogp = mlogp, motif = tmotif))
        
        
        W_g <- wilcox.test(GC_1$GC_exon, GC_2$GC_exon, alternative = "greater")
        W_l <- wilcox.test(GC_1$GC_exon, GC_2$GC_exon, alternative = "less")
        mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        # mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], type = "GC_exon", mlogp = mlogp, motif = tmotif))
        
        W_g <- wilcox.test(GC_1$GC_intron_3prime, GC_2$GC_intron_3prime, alternative = "greater")
        W_l <- wilcox.test(GC_1$GC_intron_3prime, GC_2$GC_intron_3prime, alternative = "less")
        mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        # mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], type = "GC_intron_3prime", mlogp = mlogp, motif = tmotif))
        
      } else {
        
        W_g <- wilcox.test(GC_1$Len_intron_5prime, GC_2$Len_intron_5prime, alternative = "greater")
        W_l <- wilcox.test(GC_1$Len_intron_5prime, GC_2$Len_intron_5prime, alternative = "less")
        mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        # mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], type = "Len_intron_5prime", mlogp = mlogp, motif = tmotif))
      
      
        W_g <- wilcox.test(GC_1$Len_exon, GC_2$Len_exon, alternative = "greater")
        W_l <- wilcox.test(GC_1$Len_exon, GC_2$Len_exon, alternative = "less")
        mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        # mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], type = "Len_exon", mlogp = mlogp, motif = tmotif))
      
      
        W_g <- wilcox.test(GC_1$Len_intron_3prime, GC_2$Len_intron_3prime, alternative = "greater")
        W_l <- wilcox.test(GC_1$Len_intron_3prime, GC_2$Len_intron_3prime, alternative = "less")
        mlogp <- ifelse(W_g$p.value < W_l$p.value, -log10(W_g$p.value), log10(W_l$p.value))
        # mlogp <- -log10(W_g$p.value) + log10(W_l$p.value)
        Ps <- rbind(Ps, data.frame(x = ind2sptype[i1], y = ind2sptype[i2], type = "Len_intron_3prime", mlogp = mlogp, motif = tmotif))
      
      }
    }
  }

  Ps$mlogp <- pmin(Ps$mlogp, 10)
  Ps$mlogp <- pmax(Ps$mlogp, -10)
  # Ps[Ps$x == Ps$y,"mlogp"] <- 0
  
  if (is_GC == TRUE) {
    Ps$type <- factor(Ps$type,
                      levels = c("GC_intron_5prime", "GC_exon", "GC_intron_3prime"),
                      labels = c("5' intron", "Exon", "3' intron"))
    if (is_donor == TRUE) {
      ttitle <- "GC contents, donor"
    } else {
      ttitle <- "GC contents, acceptor"
    }
  } else {
    Ps$type <- factor(Ps$type,
                      levels = c("Len_intron_5prime", "Len_exon", "Len_intron_3prime"),
                      labels = c("5' intron", "Exon", "3' intron"))   
    if (is_donor == TRUE) {
      ttitle <- "Length, donor"
    } else {
      ttitle <- "Length, acceptor"
    }
    
  }
  
  p <- ggplot(Ps, aes(x = x, y = y, fill = mlogp)) + 
    geom_tile(colour = "grey30") + 
    ggtitle(ttitle) +
    coord_flip() +
    facet_grid(.~type) +
    my_theme() +
    scale_fill_gradient2(low = "#2166ac", mid = "#ffffff", high = "#b2182b") +
    theme(axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
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
  
  # if (is_donor == TRUE) {
  #   p <- p + theme(axis.text.x = element_blank())
  # } else {
  #   p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1),
  #                  strip.background = element_blank(),
  #                  strip.text.x = element_blank(),
  #                  strip.text.y = element_blank()
  #                  )
  # }
  
    
  # labs(x = "Splicing pattern 1", y = "Splicing pattern 2", fill = "Differential score for splicing pattern 1 v.s. 2")
  
  if (is_dummy == FALSE) {
    p <- p + guides(fill = FALSE) +
      labs(x = "Splicing pattern 1", y = "Splicing pattern 2")
  } else {
    p <- p + labs(fill = "Differential score for splicing pattern 1 v.s. 2")
  }
  
  return(p)
}

##########
# for legend

# p_legend <- g_legend(get_print_info(GC, FALSE, FALSE, TRUE))

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
###########

p_GC_donor <- get_print_info(GC, TRUE, TRUE, FALSE)
p_GC_acceptor <- get_print_info(GC, FALSE, TRUE, FALSE)

# plot_grid(
#   plot_grid(p_GC_donor, p_GC_acceptor, align = "v", ncol = 2),
#   plot_grid(g_legend(p_legend_less), g_legend(p_legend_greater), nrow = 1),
#   ncol = 1, rel_heights = c(1, 0.2))
plot_grid(p_GC_donor, p_GC_acceptor, align = "v", ncol = 2)

ggsave("../figure/GC_comp.tiff", width = 18, height = 5.5, dpi = 600, units = "cm")

p_Len_donor <- get_print_info(GC, TRUE, FALSE, FALSE)
p_Len_acceptor <- get_print_info(GC, FALSE, FALSE, FALSE)

# plot_grid(
#   plot_grid(p_Len_donor, p_Len_acceptor, align = "v", ncol = 2),
#   plot_grid(g_legend(p_legend_less), g_legend(p_legend_greater), nrow = 1),
#   ncol = 1, rel_heights = c(1, 0.2))
plot_grid(p_Len_donor, p_Len_acceptor, align = "v", ncol = 2)

ggsave("../figure/Len_comp.tiff", width = 18, height = 5.5, dpi = 600, units = "cm")
  



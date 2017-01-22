library(ggplot2)
library(dplyr)

source("subscript_matome/plot_config.R")

splicing_margin <- 0.1
splicing_wide_margin <- 300

refGene <- read.table("/home/yshira/.local/lib/python2.7/site-packages/annot_utils/data/hg19/refGene.txt.gz",
                      sep = "\t", header = FALSE, stringsAsFactors = FALSE)


splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)
splicing_mutation[splicing_mutation$Splicing_Class == "intronic-alternative-5'-splice-site", "Splicing_Class"] <- "alternative-5'-splice-site"
splicing_mutation[splicing_mutation$Splicing_Class == "intronic-alternative-3'-splice-site", "Splicing_Class"] <- "alternative-3'-splice-site"
splicing_mutation[splicing_mutation$Splicing_Class == "opposite-side-intron-retention", "Splicing_Class"] <- "intron-retention"

splicing_mutation$Splicing_Class <-
  factor(splicing_mutation$Splicing_Class,
         levels = c("exon-skip", "alternative-5'-splice-site",
                    "alternative-3'-splice-site",
                    "intron-retention"),
         labels = c("Exon skip", "Alternative 5' splice site",
                    "Alternative 3' splice site", "Intron retention"))

##########
# function

get_print_info <- function(gene_symbol, mutation_key) {
  
  
  ##########
  # get gene info
  refGene_part <- refGene %>% filter(V13 == gene_symbol)
  start_target <- Inf
  end_target <- -Inf
  gene_target <- c()
  for (i in 1:nrow(refGene_part)) {
    chr_target <- refGene_part[i, 3]
    # start_target <- min(start_target, refGene_part[i, 5])
    # end_target <- max(end_target, refGene_part[i, 6])
    gene_target <- c(gene_target, refGene_part[i, 2])
  }
  ##########
  
  ##########
  # splice info
  splicing_line <- c()
  sp_mut_filt <- splicing_mutation %>% filter(Mutation_Key == mutation_key)
  
  current_y <- 1.5
  for (i in 1:nrow(sp_mut_filt)) {
    sp_start_end <- strsplit(strsplit(sp_mut_filt[i, "Splicing_Key"], ":")[[1]][2], "-")[[1]]
    sp_start <- as.numeric(sp_start_end[1])
    sp_end <- as.numeric(sp_start_end[2]) 

    if (sp_mut_filt[i, "Splicing_Class"] == "Intron retention") {
      sp_start <- sp_start - 10
      sp_end <- sp_end + 10
    }  
    start_target <- min(start_target, sp_start - splicing_wide_margin)
    end_target <- max(end_target, sp_end + splicing_wide_margin)
    
    # splicing start segment
    if (sp_start >= start_target & sp_start <= end_target) {
      tx <- sp_start
      txend <- sp_start
      ty <- current_y
      tyend <- current_y + splicing_margin
      splicing_line <- rbind(splicing_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend, splicing_class = sp_mut_filt[i, "Splicing_Class"]))
    }
    # splicing start-end segment
    if (sp_start < end_target & sp_end >=start_target) {
      tx <- max(sp_start, start_target)
      txend <- min(sp_end, end_target)
      ty <- current_y + splicing_margin
      tyend <- current_y + splicing_margin
      splicing_line <- rbind(splicing_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend, splicing_class = sp_mut_filt[i, "Splicing_Class"]))
    }
    # splicing end segment
    if (sp_end >= start_target & sp_end <= end_target) {
      tx <- sp_end
      txend <- sp_end
      ty <- current_y
      tyend <- current_y + splicing_margin
      splicing_line <- rbind(splicing_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend, splicing_class = sp_mut_filt[i, "Splicing_Class"]))
    }
    
    current_y <- current_y + 1
  }
  
  
  
  

  ##########
  # get gene coordinates 
  target_gene_info <- refGene %>% filter(V2 %in% gene_target)
  
  exon_box <- data.frame(xmin = c(), xmax = c(), ymin = c(), ymax = c())
  intron_line <- data.frame(x = c(), y = c(), xend = c(), yend = c())
  exon_intron_junction <- c()
  
  for (i in 1:nrow(target_gene_info)) {
    
    starts <- as.numeric(strsplit(target_gene_info[i,10], split=",")[[1]])
    ends <- as.numeric(strsplit(target_gene_info[i,11], split=",")[[1]])
    
    for (j in 1:length(starts)) {
      
      if (starts[j] <= end_target & ends[j] >= start_target) {
        
        # completely outside coding region
        if (ends[j] <= target_gene_info[i,7] | starts[j] >= target_gene_info[i,8]) {
          txmin <- max(starts[j], start_target)
          txmax <- min(ends[j], end_target)
          tymin <- 0.4
          tymax <- 0.6
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
        } else if (starts[j] >= target_gene_info[i,7] & ends[j] <= target_gene_info[i,8]) { # completely inside coding region
          txmin <- max(starts[j], start_target)
          txmax <- min(ends[j], end_target)
          tymin <- 0.2
          tymax <- 0.8
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
        } else {
          if (starts[j] <= target_gene_info[i,7]) {
            txmin <- max(starts[j], start_target)
            txmax <- min(target_gene_info[i,7], end_target)
            tymin <- 0.4
            tymax <- 0.6
            exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
            
            txmin <- max(target_gene_info[i,7], start_target)
            txmax <- min(ends[j], end_target)
            tymin <- 0.2
            tymax <- 0.8
            exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
          } else if (ends[j] >= target_gene_info[i,8]) {
            txmin <- max(starts[j], start_target)
            txmax <- min(target_gene_info[i,8], end_target)
            tymin <- 0.2
            tymax <- 0.8
            exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
            
            txmin <- max(target_gene_info[i,8], start_target)
            txmax <- min(ends[j], end_target)
            tymin <- 0.4
            tymax <- 0.6
            exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))         
            
          }
        }
        
      } 
      
      if (starts[j] >= start_target & starts[j] <= end_target) {
        exon_intron_junction <- c(exon_intron_junction, starts[j])
      }
      if (ends[j] >= start_target & ends[j] <= end_target) {
        exon_intron_junction <- c(exon_intron_junction, ends[j])
      }
      
      # intron line
      if (j >= 2) {
        if (ends[j - 1] <= end_target & starts[j] >= start_target) {
          if (target_gene_info[1,4] == "+") {
            tx <- max(ends[j - 1], start_target)
            txend <- min(starts[j], end_target)
          } else {
            txend <- max(ends[j - 1], start_target)
            tx <- min(starts[j], end_target)          
          }
          ty <- 0.5
          tyend <- 0.5
          intron_line <- rbind(intron_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend)) 
        }
      }
      
    }
    
  }
  
  ##########
  # get mut info
  mut_data <- strsplit(mutation_key, ',')[[1]]
  mut_line <- data.frame(x = as.numeric(mut_data[2]), y = 0.5)
  ##########
  
  print_info <- ggplot()
  if (nrow(exon_box) > 0) {
    print_info <- print_info + geom_rect(data = exon_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey80")
  }
  
  if (nrow(intron_line) > 0) {
    print_info <- print_info + geom_segment(data = intron_line, aes(x = x, xend = xend, y = y, yend = yend), colour = "gray80", arrow = arrow(length = unit(0.1, "inches"))) 
  }
  
  print_info <- print_info + theme_bw() +
    theme(legend.position = "bottom",
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_line(colour="grey50", linetype="dotted"),
          panel.border=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    scale_x_continuous(limits = c(start_target, end_target), minor_breaks = exon_intron_junction) +
    scale_y_continuous(trans = "reverse") + 
    guides(fill=FALSE) +
    labs(colour = "")
    # ggtitle(paste(gene_symbol, mutation_key))

  
  print_info <- print_info + geom_point(data = mut_line, aes(x = x, y = y), colour = "red", shape = 4)

  print_info <- print_info + geom_segment(data = splicing_line, aes(x = x, xend = xend, y = y, yend = yend, colour = splicing_class)) +
    scale_colour_manual(values = splicing_class_colour)

           
  return(print_info)                        
  
}


multiple_effect <- splicing_mutation %>% 
  group_by(Cancer_Type, Sample_Name, Gene_Symbol, Mutation_Key) %>% 
  summarize(splice_count = n()) %>%
  arrange(desc(splice_count))

multiple_effect_count <- multiple_effect %>%
  group_by(splice_count) %>% 
  summarize(count = n()) 


multiple_effect_filt <- multiple_effect %>% filter(splice_count >= 6)


if (!file.exists("../matome/multi_splice")) {
  dir.create("../matome/multi_splice")
}

for(i in 1:nrow(multiple_effect_filt)) {
  
  gene_symbol <- as.character(multiple_effect_filt[i, "Gene_Symbol"])
  mutation_key <- as.character(multiple_effect_filt[i, "Mutation_Key"])
  
  get_print_info(gene_symbol, mutation_key)
  ggsave(paste("../matome/multi_splice/", gene_symbol, ".pdf", sep =""), width = 6, height = 3.5)
  
}




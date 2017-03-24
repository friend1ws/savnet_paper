library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RColorBrewer)

source("../../../conf/plot_config.R")

##########
# setting data
splicing_class_colour <- c(
  "Exon skip" = "#fb8072",
  "Alternative 5'-ss" = "#80b1d3",
  "Alternative 3'-ss" = "#b3de69",
  "Intron retention" = "#bebada"
)

##########

gene2ref <- read.table("../db/gene2ref.txt", sep = "\t", header = FALSE)
# gene2ref <- read.table("gene2ref.txt", sep = "\t", header = FALSE)



gene2size <- gene2ref %>% 
  select(gene = V1, size = V3) %>% 
  spread(key = gene, value = size)

domain_info <- read.table("../db/domain_info.txt", 
# domain_info <- read.table("domain_info.txt", 
                          header = TRUE, sep = "\t",  
                          comment.char = "", stringsAsFactors = FALSE)

##########
# get print_gene info

get_gene_print_info <- function(ref_gene_id, start_target, end_target, dir_target, amino_size = NULL, domain_info_gene = NULL) {
  
  
  target_gene_info <- refGene %>% filter(V2 == ref_gene_id)
  
  starts <- as.numeric(strsplit(as.character(target_gene_info[10]), split=",")[[1]])
  ends <- as.numeric(strsplit(as.character(target_gene_info[11]), split=",")[[1]])
  
  start_target <- max(start_target, as.numeric(target_gene_info[5]))
  end_target <- min(end_target, as.numeric(target_gene_info[6]))
  
  gene_domain_colour <- c("grey50")
  names(gene_domain_colour) <- c("gene")
  
  ###
  # when adding amino seq info
  if (!is.null(amino_size)) {
    
    rect_data <- data.frame(
      xmin = start_target + ((end_target - start_target) / amino_size) * c(0, domain_info_gene$Pos1),
      xmax = start_target + ((end_target - start_target) / amino_size) * c(amino_size, domain_info_gene$Pos2),
      ymin = rep(0.05, 1 + nrow(domain_info_gene)),
      ymax = rep(0.95, 1 + nrow(domain_info_gene)),
      domain = c("whole", as.character(domain_info_gene$Domain))
    )
  
    if (dir_target == "-") {
      rect_data$xmin <- end_target - rect_data$xmin + start_target
      rect_data$xmax <- end_target - rect_data$xmax + start_target
    }
  
    label_data <- data.frame(
      x = (rect_data$xmin[2:nrow(rect_data)] + rect_data$xmax[2:nrow(rect_data)]) / 2,
      y = rep(0.5, nrow(domain_info_gene)),
      label = as.character(domain_info_gene$Domain)
    )
  
    domain_colour <- c("#e0e0e0", domain_info_gene$Colour, "#8c510a", "#dfc27d")

    names(domain_colour) <- c("whole", domain_info_gene$Domain, paste("ga_seg_", 1:2, sep = ""))
  
    gene_domain_colour <- c(gene_domain_colour, domain_colour)
  }
  
  
  exon_box <- data.frame(xmin = c(), xmax = c(), ymin = c(), ymax = c())
  intron_line <- data.frame(x = c(), y = c(), xend = c(), yend = c())
  exon_intron_junction <- c()
  amino_genome_seg <- data.frame(x = c(), y = c(), xend = c(), yend = c(), col = c())
  exon_num <- data.frame(x = c(), y = c(), label = c())
  
  current_coding_pos <- 0
  amino_genome_seg_color <- 1
  amino_genome_poly_group <- 1
  for (j in 1:length(starts)) {
    
    label <- paste(ifelse(target_gene_info$V4 == "+", j, length(starts) - j + 1))
    
    # completely outside coding region
    if (ends[j] <= as.numeric(target_gene_info[7]) | starts[j] >= as.numeric(target_gene_info[8])) {
      txmin <- max(starts[j], start_target)
      txmax <- min(ends[j], end_target)
      tymin <- 2 + 0.3
      tymax <- 2 + 0.7
      if (starts[j] <= end_target & ends[j] >= start_target) {
        exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
        # exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 3.25, label = label))
      }
    # completely inside coding region
    } else if (starts[j] >= as.numeric(target_gene_info[7]) & ends[j] <= as.numeric(target_gene_info[8])) { 
      txmin <- max(starts[j], start_target)
      txmax <- min(ends[j], end_target)
      tymin <- 2 + 0.05
      tymax <- 2 + 0.95
      if (starts[j] <= end_target & ends[j] >= start_target) {
        exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
        exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 3.25, label = label))
      }
      
      if (!is.null(amino_size)) {
        cx1 <- start_target + ((end_target - start_target) / amino_size) * current_coding_pos
        current_coding_pos <- current_coding_pos + (ends[j] - starts[j]) / 3
        cx2 <- start_target + ((end_target - start_target) / amino_size) * (current_coding_pos - 1)
        if (starts[j] <= end_target & ends[j] >= start_target) {
          amino_genome_seg <- rbind(amino_genome_seg,
                                  data.frame(x = c(starts[j], ends[j], cx2, cx1), 
                                             y = c(2, 2, 1, 1),
                                             group = rep(paste("ga_seg_", amino_genome_poly_group, sep = ""), 4),
                                             col = rep(paste("ga_seg_", amino_genome_seg_color, sep = ""), 4)
                                             )
                                  )
        }
      }
      
    } else {
      if (starts[j] <= as.numeric(target_gene_info[7])) {
        txmin <- max(starts[j], start_target)
        txmax <- min(as.numeric(target_gene_info[7]), end_target)
        tymin <- 2 + 0.3
        tymax <- 2 + 0.7
        if (starts[j] <= end_target & ends[j] >= start_target) {
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
          # exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 3.25, label = label))
        }
        
        txmin <- max(as.numeric(target_gene_info[7]), start_target)
        txmax <- min(ends[j], end_target)
        tymin <- 2 + 0.05
        tymax <- 2 + 0.95
        if (starts[j] <= end_target & ends[j] >= start_target) {
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
          exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 3.25, label = label))
        }
        
        if (!is.null(amino_size)) {
          
          cx1 <- start_target + ((end_target - start_target) / amino_size) * current_coding_pos
          current_coding_pos <- current_coding_pos + (ends[j] - as.numeric(target_gene_info[7])) / 3
          cx2 <- start_target + ((end_target - start_target) / amino_size) * (current_coding_pos - 1)
          
          if (starts[j] <= end_target & ends[j] >= start_target) {
            amino_genome_seg <- rbind(amino_genome_seg,
                                      data.frame(x = c(as.numeric(target_gene_info[7]), ends[j], cx2, cx1),
                                                 y = c(2, 2, 1, 1),
                                                 group = rep(paste("ga_seg_", amino_genome_poly_group, sep = ""), 4),
                                                 col = rep(paste("ga_seg_", amino_genome_seg_color, sep = ""), 4)
                                               )
                                    )
          }
         }
        
      } else if (ends[j] >= as.numeric(target_gene_info[8])) {
        txmin <- max(starts[j], start_target)
        txmax <- min(as.numeric(target_gene_info[8]), end_target)
        tymin <- 2 + 0.05
        tymax <- 2 + 0.95
        if (starts[j] <= end_target & ends[j] >= start_target) {
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
          exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 3.25, label = label))
        }
          
        if (!is.null(amino_size)) {
            
          cx1 <- start_target + ((end_target - start_target) / amino_size) * current_coding_pos
          current_coding_pos <- current_coding_pos + (as.numeric(target_gene_info[8]) - starts[j]) / 3
          cx2 <- start_target + ((end_target - start_target) / amino_size) * (current_coding_pos - 1)
            
          if (starts[j] <= end_target & ends[j] >= start_target) {
            amino_genome_seg <- rbind(amino_genome_seg,
                                      data.frame(x = c(starts[j], as.numeric(target_gene_info[8]), cx2, cx1), 
                                                 y = c(2, 2, 1, 1), 
                                                 group = rep(paste("ga_seg_", amino_genome_poly_group, sep = ""), 4),
                                                 col = rep(paste("ga_seg_", amino_genome_seg_color, sep = ""), 4)
                                      )
            )
          }
          
        }
        
        txmin <- max(as.numeric(target_gene_info[8]), start_target)
        txmax <- min(ends[j], end_target)
        tymin <- 2 + 0.3
        tymax <- 2 + 0.7
        if (starts[j] <= end_target & ends[j] >= start_target) {
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))  
          # exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 3.25, label = label))
        }
        
      }
    
    }
    
    amino_genome_seg_color <- amino_genome_seg_color + 1
    if (amino_genome_seg_color > 2) amino_genome_seg_color <- 1
    amino_genome_poly_group <- amino_genome_poly_group + 1
    # if (amino_genome_poly_group > 8) amino_genome_poly_group <- 1

  
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
        ty <- 2 + 0.5
        tyend <- 2 + 0.5
        intron_line <- rbind(intron_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend)) 
      }
    }
    
  }
  
  
  # amino_genome_seg$xend <- start_target + ((end_target - start_target) / current_coding_pos) * amino_genome_seg$xend
  amino_genome_seg$x <- pmax(amino_genome_seg$x, start_target)
  amino_genome_seg$x <- pmin(amino_genome_seg$x, end_target)  
    
  p_gene <- ggplot()
  if (nrow(exon_box) > 0) {
    p_gene <- p_gene + geom_rect(data = exon_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill= "gene"))
    p_gene <- p_gene + geom_text(data = exon_num, aes(x = x, y = y, label = label, colour = "gene"), size = 1.8)
  }
  
  if (!is.null(amino_size)) {
    p_gene <- p_gene + geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = domain)) +
      geom_text(data = label_data, aes(x = x, y = y, label = label), size = 1.8)
  
    # p_gene <- p_gene + geom_segment(data = amino_genome_seg, 
    #                                 aes(x = x, xend = xend, y = y, yend = yend, color = col),
    #                                 linetype = "dashed", size = 0.2)
    p_gene <- p_gene + geom_polygon(data = amino_genome_seg,
                                    aes(x = x, y = y, group = group, fill = col), alpha = 0.5)
    
  }
  
  if (nrow(intron_line) > 0) {
    p_gene <- p_gene + geom_segment(data = intron_line, aes(x = x, xend = xend, y = y, yend = yend, colour = "gene"), 
                                    size = 0.3, arrow = arrow(length = unit(0.03, "inches"))) 
  }
  
  p_gene <- p_gene + theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.border=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines")) +
    scale_y_continuous(trans = "reverse") + 
    scale_fill_manual(values = gene_domain_colour) +
    scale_colour_manual(values = gene_domain_colour) +
    guides(fill=FALSE) +
    guides(colour = FALSE) + 
    guides(linetype = FALSE) + 
    guides(size = FALSE) 
  
  if (dir_target == "-") {
    p_gene <- p_gene + scale_x_continuous(limits = c(end_target, start_target), trans = "reverse", minor_breaks = exon_intron_junction)
  } else {
    p_gene <- p_gene + scale_x_continuous(limits = c(start_target, end_target), minor_breaks = exon_intron_junction) 
  } 
  
  return(list(p_gene, exon_box, intron_line, exon_intron_junction))
  
}



##########


##########
get_gsm_print_info <- function(omega_info_input, gene_symbol, ref_gene_id, count_thres, start_target, end_target, dir_target, 
                               mut_type_str = "", exon_intron_junction = c(), mut_margin = 300) {
  
  
  splicing_mut_class_colour <- c(splicing_class_colour, c("donor" = "#80cdc1", "acceptor" = "#91cf60"))
  
  
  omega_info_input <- omega_info_input %>% filter(Gene_Symbol == gene_symbol)

  intron_ind <- grepl("Intron retention", omega_info_input$Splicing_Class)
   
  omega_info_input$Splicing_Key[intron_ind] <- 
    unlist(lapply(omega_info_input$Splicing_Key[intron_ind], 
                  function(x) {get_intron_key(x, ref_gene_id)}))
  
  omega_info_input$Is_Inframe2 <- omega_info_input$Is_Inframe
  omega_info_input$Is_Inframe2[!intron_ind] <- 
    unlist(lapply(omega_info_input$Splicing_Key[!intron_ind], 
                  function(x) {check_inframe(x, ref_gene_id)}))
  
  

  motif_level_mut_count <- omega_info_input %>% 
    group_by(Motif_Pos) %>% 
    summarize(count = n()) %>%
    filter(count >= count_thres)

  # if (gene_symbol == "CDKN2A") { 
  #   print(as.data.frame(motif_level_mut_count))
  # }

  sp_mut_filt <- omega_info_input %>% 
    filter(Motif_Pos %in% motif_level_mut_count$Motif_Pos) %>%
    group_by(Gene_Symbol, Motif_Pos, Mutation_Type, Splicing_Key, Splicing_Class, Is_Inframe2) %>%
    summarize(count = n()) 
  
  if (mut_type_str != "") {
    sp_mut_filt <- sp_mut_filt %>% filter(grepl(mut_type_str, Mutation_Type))
  }
  
  
  if (dir_target == "-") {
    sp_mut_filt <- sp_mut_filt %>% arrange(desc(Motif_Pos), desc(count))
  } else {
    sp_mut_filt <- sp_mut_filt %>% arrange(Motif_Pos, desc(count))
  }
  
  # if (gene_symbol == "CDKN2A") {  
  #   print("sp_mut_filt_start")
  #   print(as.data.frame(sp_mut_filt))
  #   print("sp_mut_filt_end")
  # }

  mut_line <- c()
  mut_seg <- c()
  splicing_line <- c()
  tmp_splicing_line <- c()
  current_y_minor <- 0
  current_max_x_pos <- c()
  
  splicing_margin <- 0.5
  tmp_mutation_num <- 0
  
  if (nrow(sp_mut_filt) == 0) {
    return(list(ggplot(), splicing_line, mut_line))
  }
  
  temp_mut <- data.frame(Motif_Pos = c(""))
  temp_mut_pos <- ""
  for (i in 1:nrow(sp_mut_filt)) {
    
    mut_data <- strsplit(strsplit(sp_mut_filt$Motif_Pos[i], ":")[[1]][2], "-")[[1]]
    mut_pos <- as.numeric(mut_data[1])
    if (mut_pos < start_target | mut_pos > end_target) next
    
    if (mut_pos != temp_mut_pos) {
      
      if (temp_mut_pos != "") {

        max_x_pos <- ifelse(dir_target == "+", 
                            max(tmp_splicing_line$xend) + mut_margin,
                            min(tmp_splicing_line$x) - mut_margin)
        
        if (length(current_max_x_pos) == 0) {
          current_y <- 0
          current_max_x_pos[1] <- max_x_pos
        } else {
          prev_y_flag <- FALSE
          for (j in 1:length(current_max_x_pos)) {
            # print(c(min(tmp_splicing_line$x), current_max_x_pos[j]))
            if ((dir_target == "+" & min(tmp_splicing_line$x) > current_max_x_pos[j]) |
                (dir_target == "-" & max(tmp_splicing_line$xend) < current_max_x_pos[j])) {
              current_y <- 9 * (j - 1)
              current_max_x_pos[j] <- max_x_pos
              prev_y_flag <- TRUE
              break
            }
          }
        
          if (prev_y_flag == FALSE) {
            current_y <- 9 * j
            current_max_x_pos <- c(current_max_x_pos, max_x_pos)
          }
        }
        
        mut_line <- rbind(mut_line, 
                          data.frame(x = as.numeric(temp_mut_pos), y = current_y + 5,
                                     motif_type = ifelse(grepl("donor", temp_mut$Mutation_Type), "donor", "acceptor"),
                                     mutation_type = ifelse(grepl("disruption", temp_mut$Mutation_Type), "disruption", "creation"),
                                     mutation_num = ifelse(temp_mut_count > 20, 20, temp_mut_count),
                                     count_label = ifelse(temp_mut_count >= 2, as.character(temp_mut_count), "")))
        
        mut_seg <- rbind(mut_seg,
                         data.frame(x = as.numeric(temp_mut_pos), y = current_y + 5,
                                    xend = as.numeric(temp_mut_pos), yend = current_y))

        tmp_splicing_line$y <- tmp_splicing_line$y + current_y
        tmp_splicing_line$yend <- tmp_splicing_line$yend + current_y
        
        splicing_line <- rbind(splicing_line, tmp_splicing_line)
      }
      
      temp_mut <- sp_mut_filt[i,]
      temp_mut_pos <- mut_pos
      tmp_splicing_line <- c()
      current_y_minor <- 0
      temp_mut_count <- 0
      
    }
    
    
    # get mut info
    
    temp_mut_count <- temp_mut_count + as.numeric(sp_mut_filt[i,"count"])
    sp_start_end <- strsplit(strsplit(as.character(sp_mut_filt[i, "Splicing_Key"]), ":")[[1]][2], "-")[[1]]
    sp_start <- as.numeric(sp_start_end[1])
    sp_end <- as.numeric(sp_start_end[2]) 

    
    # splicing start segment
    if (sp_start >= start_target & sp_start <= end_target) {
      tx <- sp_start
      txend <- sp_start
      ty <- current_y_minor
      tyend <- current_y_minor + splicing_margin
      tmp_splicing_line <- rbind(tmp_splicing_line, 
                             data.frame(x = tx, xend = txend, 
                                        y = ty, yend = tyend, 
                                        splicing_class = sp_mut_filt$Splicing_Class[i], 
                                        is_inframe = ifelse(sp_mut_filt$Is_Inframe2[i] == "in-frame", "In frame", "Frameshift")))
    }
    
    # splicing start-end segment
    if (sp_start < end_target & sp_end >=start_target) {
      tx <- max(sp_start, start_target)
      txend <- min(sp_end, end_target)
      ty <- current_y_minor + splicing_margin
      tyend <- current_y_minor + splicing_margin
      tmp_splicing_line <- rbind(tmp_splicing_line, 
                             data.frame(x = tx, xend = txend, 
                                        y = ty, yend = tyend,
                                        splicing_class = sp_mut_filt$Splicing_Class[i], 
                                        is_inframe = ifelse(sp_mut_filt$Is_Inframe2[i] == "in-frame", "In frame", "Frameshift")))
    }
    
    # splicing end segment
    if (sp_end >= start_target & sp_end <= end_target) {
      tx <- sp_end
      txend <- sp_end
      ty <- current_y_minor
      tyend <- current_y_minor + splicing_margin
      tmp_splicing_line <- rbind(tmp_splicing_line, 
                             data.frame(x = tx, xend = txend,
                                        y = ty, yend = tyend, 
                                        splicing_class = sp_mut_filt$Splicing_Class[i], 
                                        is_inframe = ifelse(sp_mut_filt$Is_Inframe2[i] == "in-frame", "In frame", "Frameshift")))
    }
    
    current_y_minor <- current_y_minor + 1 
   
    
  }

  
  # last processing
  if (temp_mut_pos != "") {
    
    max_x_pos <- ifelse(dir_target == "+", 
                        max(tmp_splicing_line$xend) + mut_margin,
                        min(tmp_splicing_line$x) - mut_margin)
    
    if (length(current_max_x_pos) == 0) {
      current_y <- 0
      current_max_x_pos[1] <- max_x_pos
    } else {
      prev_y_flag <- FALSE
      for (j in 1:length(current_max_x_pos)) {
        if ((dir_target == "+" & min(tmp_splicing_line$x) > current_max_x_pos[j]) |
            (dir_target == "-" & max(tmp_splicing_line$xend) < current_max_x_pos[j])) {
          current_y <- 9 * (j - 1)
          current_max_x_pos[j] <- max_x_pos
          prev_y_flag <- TRUE
          break
        }
      }
      
      if (prev_y_flag == FALSE) {
        current_y <- 9 * j
        current_max_x_pos <- c(current_max_x_pos, max_x_pos)
      }
    }
    
    mut_line <- rbind(mut_line, 
                      data.frame(x = as.numeric(temp_mut_pos), y = current_y + 5,
                                 motif_type = ifelse(grepl("donor", temp_mut$Mutation_Type), "donor", "acceptor"),
                                 mutation_type = ifelse(grepl("disruption", temp_mut$Mutation_Type), "disruption", "creation"),
                                 mutation_num = ifelse(temp_mut$count > 20, 20, temp_mut$count),
                                 count_label = ifelse(temp_mut$count >= 2, as.character(temp_mut$count), "")))
    
    mut_seg <- rbind(mut_seg,
                     data.frame(x = as.numeric(temp_mut_pos), y = current_y + 5,
                                xend = as.numeric(temp_mut_pos), yend = current_y))
    
    tmp_splicing_line$y <- tmp_splicing_line$y + current_y
    tmp_splicing_line$yend <- tmp_splicing_line$yend + current_y
    
    splicing_line <- rbind(splicing_line, tmp_splicing_line)
  }
  
 
 
  p_gsm <- ggplot() + theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_line(colour="grey60", linetype="dotted", size = 0.25),
          panel.border=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines")) +
    scale_y_continuous(trans = "reverse", limits = c(max(mut_line$y) + 1.5, 0)) + 
    guides(fill=FALSE) +
    guides(colour = FALSE) + 
    guides(linetype = FALSE) + 
    guides(size = FALSE) +
    guides(shape = FALSE)
  
  if (dir_target == "-") {
    p_gsm <- p_gsm + scale_x_continuous(limits = c(end_target, start_target), trans = "reverse", minor_breaks = exon_intron_junction) 
  } else {
    p_gsm <- p_gsm + scale_x_continuous(limits = c(start_target, end_target), minor_breaks = exon_intron_junction) 
  } 
  
  p_gsm <- p_gsm + geom_segment(data = mut_seg,
                                aes(x = x, y = y, xend = xend, yend = yend), colour = "grey50", size = 0.5)
  
  p_gsm <- p_gsm + geom_point(data = mut_line, 
                              aes(x = x, y = y, shape = mutation_type, fill = motif_type, size = mutation_num, colour = motif_type), alpha = 0.95) +
    scale_shape_manual(values = c(disruption = 21, creation = 24)) +
    scale_size(range = c(1, 3), limits = c(1, 20))
  
  p_gsm <- p_gsm + geom_text(data = mut_line, 
                             aes(x = x, y = y, label = count_label), colour = "gray20", size = 1.8) 
  
  
  p_gsm <- p_gsm + geom_segment(data = splicing_line, 
                                aes(x = x, xend = xend, y = y, yend = yend, 
                                    colour = splicing_class, linetype = is_inframe), size = 0.5) +
    scale_colour_manual(values = splicing_mut_class_colour) +
    scale_fill_manual(values = splicing_mut_class_colour) +
    scale_linetype_manual(values = c("In frame" = "solid", "Frameshift" = "longdash"))
  
  return(list(p_gsm, splicing_line, mut_line))
  
}

##########
# get intron pos for opposite intron retention
get_intron_key <- function(motif_pos, ref_gene_id) {

  motif_pos_sp1 <- strsplit(motif_pos, ":")
  motif_pos_sp2 <- strsplit(unlist(lapply(motif_pos_sp1, '[', 2)), ",")
  motif_pos_sp3 <- strsplit(unlist(lapply(motif_pos_sp2, '[', 1)), "-")

  tchr <- unlist(lapply(motif_pos_sp1, '[', 1))
  tstart <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 1)))
  tend <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 2)))
 
  target_gene_info <- refGene %>% filter(V2 == ref_gene_id)
  exon_starts <- as.numeric(strsplit(as.character(target_gene_info[10]), split=",")[[1]])
  exon_ends <- as.numeric(strsplit(as.character(target_gene_info[11]), split=",")[[1]])
  
  exon_ind1 <- which(abs(exon_starts - tstart) < 5)
  exon_ind2 <- which(abs(exon_ends - tstart) < 5)
  
  if (length(exon_ind1) > 0) {
    return(paste(tchr, paste(exon_ends[exon_ind1 - 1], exon_starts[exon_ind1], sep = "-"), sep = ":"))
  } else if (length(exon_ind2) > 0) {
    return(paste(tchr, paste(exon_ends[exon_ind2], exon_starts[exon_ind2 + 1], sep = "-"), sep = ":"))
  } else {  
    return(motif_pos)
  }
}

##########
# check_inframe
check_inframe <- function(splicing_key, ref_gene_id) {
  
  splicing_key_sp1 <- strsplit(splicing_key, ":")
  splicing_key_sp2 <- strsplit(unlist(lapply(splicing_key_sp1, '[', 2)), ",")
  splicing_key_sp3 <- strsplit(unlist(lapply(splicing_key_sp2, '[', 1)), "-")
  
  tchr <- unlist(lapply(splicing_key_sp1, '[', 1))
  tstart <- as.numeric(unlist(lapply(splicing_key_sp3, '[', 1)))
  tend <- as.numeric(unlist(lapply(splicing_key_sp3, '[', 2)))
  
  target_gene_info <- refGene %>% filter(V2 == ref_gene_id)
  exon_starts <- as.numeric(strsplit(as.character(target_gene_info[10]), split=",")[[1]]) + 1
  exon_ends <- as.numeric(strsplit(as.character(target_gene_info[11]), split=",")[[1]])
  
  coding_starts <- pmax(exon_starts, target_gene_info$V7)
  coding_ends <- pmin(exon_ends, target_gene_info$V8)
  
  spliced_coding_size <- 0
  smargin <- 30
  for (i in 1:length(coding_starts)) {
    # splice start site overlap with coding region
    if (tstart >= coding_starts[i] - smargin & tstart <= coding_ends[i] + smargin) {
      # splice end site overlap with coding region
      if (tend >= coding_starts[i] - smargin & tend <= coding_ends[i] + smargin) {
        spliced_coding_size <- spliced_coding_size + tend - tstart + 1
      } else {
        spliced_coding_size <- spliced_coding_size + coding_ends[i] - tstart + 1
      }
    # splice end site overlap with coding region
    } else if (tend >= coding_starts[i] - smargin & tend <= coding_ends[i] + smargin) {
      spliced_coding_size <- spliced_coding_size + tend - coding_starts[i] + 1
    } else if (coding_starts[i] >= tstart & coding_ends[i] <= tend) {
      # spliced region include coding region
        spliced_coding_size <- spliced_coding_size + coding_ends[i] - coding_starts[i] + 1
    } 
    
  }
   
  return(ifelse(spliced_coding_size %% 3 == 0 & spliced_coding_size != 0, "in-frame", "---"))
}

##########

##########
# processing necessay data

refGene <- read.table("../../../db/refGene/refGene.txt.gz", header = FALSE, stringsAsFactors = FALSE)
# refGene <- read.table("refGene.txt.gz", header = FALSE, stringsAsFactors = FALSE)


omega_info <- read.table("../../output/rescue/TCGA.savnet.with_rescued.result.txt", 
# omega_info <- read.table("TCGA.savnet.with_rescued.result.txt", 
                         sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE)


omega_info[omega_info$Splicing_Class == "intronic-alternative-5'-splice-site", "Splicing_Class"] <- "alternative-5'-splice-site"
omega_info[omega_info$Splicing_Class == "intronic-alternative-3'-splice-site", "Splicing_Class"] <- "alternative-3'-splice-site"
omega_info[omega_info$Splicing_Class == "opposite-side-intron-retention", "Splicing_Class"] <- "intron-retention"

omega_info$Splicing_Class <-
  factor(omega_info$Splicing_Class,
         levels = c("exon-skip", "alternative-5'-splice-site",
                    "alternative-3'-splice-site", "intron-retention"),
         labels = c("Exon skip", "Alternative 5'-ss",
                    "Alternative 3'-ss", "Intron retention"))

omega_info_proc <- omega_info %>% group_by(Sample_Name, Mutation_Key) %>% 
  summarize(splicing_key_paste = paste(Splicing_Key, collapse = ";"), 
            read_num_paste = paste(Supporting_Read_Num, collapse = ","))


uniq_key <- unlist(
  lapply(
    1:nrow(omega_info_proc),
    function(x) {
      treads <- as.numeric(strsplit(omega_info_proc$read_num_paste[x], ",")[[1]])
      tind <- which(treads == max(treads))[1]
      
      return(paste(omega_info_proc$Sample_Name[x], 
                   omega_info_proc$Mutation_Key[x], 
                   strsplit(omega_info_proc$splicing_key_paste[x], ";")[[1]][tind]))
    }
  )
)

omega_info_uniq <- omega_info  %>% 
  mutate(check_key = paste(Sample_Name, Mutation_Key, Splicing_Key)) %>% 
  filter(check_key %in% uniq_key) 

##########



print_prof <- function(gene_symbol, ref_gene_id, start_target, end_target, dir_target, mut_margin = 300) {
  
  if (gene_symbol != "CDKN2A") {
    gene_print_info <- get_gene_print_info(ref_gene_id, start_target, end_target, dir_target, 
                                           as.numeric(gene2size[gene_symbol]), domain_info %>% filter(Gene_Symbol == gene_symbol))
  } else {
    gene_print_info <- gene_print_info_CDKN2A()
  }
  
  
  
  p_gene <- gene_print_info[[1]]
  exon_intron_junction <- gene_print_info[[4]]

  gsm_print_info_d <- get_gsm_print_info(omega_info_uniq, gene_symbol, ref_gene_id, 1, start_target, end_target, dir_target, "disruption", exon_intron_junction, mut_margin)
  gsm_print_info_c <- get_gsm_print_info(omega_info_uniq, gene_symbol, ref_gene_id, 1, start_target, end_target, dir_target, "creation", exon_intron_junction, mut_margin)

  p_gsm_d <- gsm_print_info_d[[1]]
  p_gsm_c <- gsm_print_info_c[[1]]

  if (gene_symbol != "CDKN2A") {
    xtitle <- ggdraw() + draw_label(paste(gene_symbol, paste("(", ref_gene_id, ", ", as.numeric(gene2size[gene_symbol]), "aa", ")", sep = "")), size = 7)
  } else {
    xtitle <- ggdraw() + draw_label("CDKN2A(NM_000077,NM_058195, 156aa,132aa)", size = 7)
  }
 
  ylabel_dummy <- ggdraw() + draw_label("", angle = 90) 
  ylabel_d <- ggdraw() + draw_label("Disruption", angle = 90)
  ylabel_c <- ggdraw() + draw_label("Creation", angle = 90)
 

  if ( !is.null(gsm_print_info_d[[3]]) & !is.null(gsm_print_info_c[[3]]) ) {
  
    mut_y_d <- max(gsm_print_info_d[[3]]$y)
    mut_y_c <- max(gsm_print_info_c[[3]]$y)
    if (mut_y_d < 10 | mut_y_c < 10)  {
      ylabel_d <- ylabel_dummy
      ylabel_c <- ylabel_dummy
    }
 
    plot_grid(xtitle, plot_grid(ylabel_dummy, p_gene, ncol = 2, rel_widths = c(0.04, 0.96)),
              plot_grid(ylabel_d, p_gsm_d, ncol = 2, rel_widths = c(0.04, 0.96)),
              plot_grid(ylabel_c, p_gsm_c, rel_widths = c(0.04, 0.96)),
              ncol = 1, align = "v", rel_heights = c(1, 3, 0.15 * c(mut_y_d, mut_y_c) + 0.6))
  
  
  } else if ( !is.null(gsm_print_info_d[[3]]) ) {
  
    mut_y_d <- max(gsm_print_info_d[[3]]$y)
    if (mut_y_d < 10) ylabel_d <- ylabel_dummy

    plot_grid(xtitle, plot_grid(ylabel_dummy, p_gene, ncol = 2, rel_widths = c(0.04, 0.96)),
              plot_grid(ylabel_d, p_gsm_d, ncol = 2, rel_widths = c(0.04, 0.96)),
              ncol = 1, align = "v", rel_heights = c(1, 3, 0.15 * c(mut_y_d) + 0.6))
  
  } else if ( !is.null(gsm_print_info_c[[3]]) ) {
  
    mut_y_c <- max(gsm_print_info_c[[3]]$y)
    if (mut_y_c < 10) ylabel_c <- ylabel_dummy
    
    plot_grid(plot_grid(ylabel_dummy, p_gene, ncol = 2, rel_widths = c(0.04, 0.96)),
              plot_grid(ylabel_c, p_gsm_c, rel_widths = c(0.04, 0.96)),
              ncol = 1, align = "v", rel_heights = c(1, 3, 0.15 * c(mut_y_c) + 0.6))
  
  
  } else {
  
    p_gene
  
  }
  
  ggsave(paste("../figure/", gene_symbol, "_prof.tiff", sep = ""), width = 10, height = 0.6 + 2.4 + 0.08 *(mut_y_d + mut_y_c), dpi = 600, units = "cm")
  

}


gene_print_info_CDKN2A <- function() {
  
  gene_domain_colour <- c("grey50")
  names(gene_domain_colour) <- c("gene")
  
  gene_print_info_CDKN2A_1 <- get_gene_print_info("NM_000077", 21968227 - 300, 21994330 + 300, "-", NULL, NULL)
  gene_print_info_CDKN2A_2 <- get_gene_print_info("NM_058195", 21968227 - 300, 21994330 + 300, "-", NULL, NULL)
  
  gene_print_info_CDKN2A_1[[2]]$ymin <- gene_print_info_CDKN2A_1[[2]]$ymin - 1
  gene_print_info_CDKN2A_1[[2]]$ymax <- gene_print_info_CDKN2A_1[[2]]$ymax - 1
  
  gene_print_info_CDKN2A_1[[3]]$y <- gene_print_info_CDKN2A_1[[3]]$y - 1
  gene_print_info_CDKN2A_1[[3]]$yend <- gene_print_info_CDKN2A_1[[3]]$yend - 1
  
  exon_box <- rbind(gene_print_info_CDKN2A_1[[2]], gene_print_info_CDKN2A_2[[2]])
  intron_line <- rbind(gene_print_info_CDKN2A_1[[3]], gene_print_info_CDKN2A_2[[3]])
  exon_intron_junction <- unique(c(gene_print_info_CDKN2A_1[[4]], gene_print_info_CDKN2A_2[[4]]))
  
  p_gene <- ggplot()
  p_gene <- p_gene + geom_rect(data = exon_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill= "gene"))
  # p_gene <- p_gene + geom_text(data = exon_num, aes(x = x, y = y, label = label, colour = "gene"), size = 3)

  p_gene <- p_gene + geom_segment(data = intron_line, aes(x = x, xend = xend, y = y, yend = yend, colour = "gene"), size = 0.5, arrow = arrow(length = unit(0.02, "inches"))) 

  
  p_gene <- p_gene + theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.border=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines")) +
    scale_y_continuous(trans = "reverse") + 
    scale_fill_manual(values = gene_domain_colour) +
    scale_colour_manual(values = gene_domain_colour) +
    guides(fill=FALSE) +
    guides(colour = FALSE) + 
    guides(linetype = FALSE) + 
    guides(size = FALSE) 
  
    p_gene <- p_gene + scale_x_continuous(limits = c(21994330 + 300, 21968227 - 300), trans = "reverse", minor_breaks = exon_intron_junction)
  
    return(list(p_gene, exon_box, intron_line, exon_intron_junction))
}


##########

print_prof("TP53", "NM_000546", 7572926 - 300, 7579912 + 300, "-", 50)

print_prof("PIK3R1", "NM_181523", 67588086 - 300, 67593429 + 300, "+")

print_prof("CDKN2A", "NM_000077", 21968227 - 300, 21994330 + 300, "-")

print_prof("GATA3", "NM_002051", 8105958  - 300, 8115986 + 300, "+")





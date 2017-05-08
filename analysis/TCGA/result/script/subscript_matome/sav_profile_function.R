library(ggplot2)
library(dplyr)
library(cowplot)

source("../../../conf/plot_config.R")

splicing_margin <- 0.2
splicing_wide_margin <- 300

refGene <- read.table("../../../db/refGene/refGene.txt.gz",
                      sep = "\t", header = FALSE, stringsAsFactors = FALSE)

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

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
# function

get_print_info <- function(splicing_mutation, sample_name, gene_symbol, mutation_key, major_sp = "") {


  start_target <- Inf
  end_target <- -Inf
  
  # get gene info
  tmp_ref_gene <- refGene %>% filter(V13 == gene_symbol)
  gene_ind <- order(as.numeric(sub("NM_", "", tmp_ref_gene$V2)))[1]
  target_gene_info <- tmp_ref_gene[gene_ind,]
  
  ##########
  # splice info
  sp_mut_filt <- splicing_mutation %>% filter(Mutation_Key == mutation_key, Sample_Name == sample_name)
  
  intron_ind <- grepl("Intron retention", sp_mut_filt$Splicing_Class)
  
  sp_mut_filt$Splicing_Key[intron_ind] <- 
    unlist(lapply(sp_mut_filt$Splicing_Key[intron_ind], 
                  function(x) {get_intron_key(x, target_gene_info$V2)}))
  
  sp_mut_filt <- sp_mut_filt %>% 
    select(Cancer_Type, Gene_Symbol, Sample_Name, Mutation_Key, 
           Motif_Pos, Mutation_Type, Splicing_Key, Splicing_Class, Is_Inframe, Supporting_Read_Num) %>% 
    group_by(Cancer_Type, Gene_Symbol, Sample_Name, Mutation_Key,
             Motif_Pos, Mutation_Type, Splicing_Key, Splicing_Class, Is_Inframe) %>%
    summarize(mSupporting_Read_Num = max(Supporting_Read_Num)) %>%
    arrange(desc(mSupporting_Read_Num)) 

   sp_mut_filt <- data.frame(sp_mut_filt)
  print(sp_mut_filt)


  # determine the start and end target position  
  for (i in 1:nrow(sp_mut_filt)) {
    sp_start_end <- strsplit(strsplit(as.character(sp_mut_filt[i, "Splicing_Key"]), ":")[[1]][2], "-")[[1]]
    sp_start <- as.numeric(sp_start_end[1])
    sp_end <- as.numeric(sp_start_end[2])

    # start_target <- min(start_target, sp_start - splicing_wide_margin)
    # end_target <- max(end_target, sp_end + splicing_wide_margin)
    start_target <- min(start_target, sp_start)
    end_target <- max(end_target, sp_end)
  }

  target_size_raw <- (end_target - start_target)

  start_target <- start_target - max(target_size_raw / 10, splicing_wide_margin)
  end_target <- end_target + max(target_size_raw / 10, splicing_wide_margin)
  sr_num_margin <- target_size_raw / 30
  

  current_y <- 2.5
  splicing_line <- c()
  sr_num <- data.frame(x = c(), y = c(), label = c())
  major_simbol <- data.frame(x = c(), y = c())

  for (i in 1:nrow(sp_mut_filt)) {
 
    sp_start_end <- strsplit(strsplit(as.character(sp_mut_filt[i, "Splicing_Key"]), ":")[[1]][2], "-")[[1]]
    sp_start <- as.numeric(sp_start_end[1])
    sp_end <- as.numeric(sp_start_end[2])
  
    # splicing start segment
    if (sp_start >= start_target & sp_start <= end_target) {
      tx <- sp_start
      txend <- sp_start
      ty <- current_y
      tyend <- current_y + splicing_margin
      splicing_line <- rbind(splicing_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend, 
                                                       splicing_class = as.character(sp_mut_filt[i, "Splicing_Class"]),
                                                       is_inframe = as.character(sp_mut_filt[i, "Is_Inframe"])))
      if (target_gene_info$V4 == "-") {
        sr_num <- rbind(sr_num, data.frame(x = tx - sr_num_margin, y = current_y + splicing_margin, label = as.character(sp_mut_filt[i, "mSupporting_Read_Num"])))
      }
      if (as.character(sp_mut_filt[i, "Splicing_Key"]) == major_sp & target_gene_info$V4 == "+") {
        major_simbol <- data.frame(x = tx - sr_num_margin, y = current_y + splicing_margin)
      }
    }
    # splicing start-end segment
    if (sp_start < end_target & sp_end >=start_target) {
      tx <- max(sp_start, start_target)
      txend <- min(sp_end, end_target)
      ty <- current_y + splicing_margin
      tyend <- current_y + splicing_margin
      splicing_line <- rbind(splicing_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend, 
                                                       splicing_class = as.character(sp_mut_filt[i, "Splicing_Class"]),
                                                       is_inframe = as.character(sp_mut_filt[i, "Is_Inframe"])))
    }
    # splicing end segment
    if (sp_end >= start_target & sp_end <= end_target) {
      tx <- sp_end
      txend <- sp_end
      ty <- current_y
      tyend <- current_y + splicing_margin
      splicing_line <- rbind(splicing_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend, 
                                                       splicing_class = as.character(sp_mut_filt[i, "Splicing_Class"]),
                                                       is_inframe = as.character(sp_mut_filt[i, "Is_Inframe"])))
      if (target_gene_info$V4 == "+") {
        sr_num <- rbind(sr_num, data.frame(x = txend + sr_num_margin, y = current_y + splicing_margin, label = as.character(sp_mut_filt[i, "mSupporting_Read_Num"])))
      }
      if (as.character(sp_mut_filt[i, "Splicing_Key"]) == major_sp & target_gene_info$V4 == "-") {
        major_simbol <- data.frame(x = tx + sr_num_margin, y = current_y + splicing_margin)
      }

    }
    
    current_y <- current_y + 1
  }
  

  ##########
  # get gene coordinates 
  

  exon_num <- data.frame(x = c(), y = c(), label = c())
  exon_box <- data.frame(xmin = c(), xmax = c(), ymin = c(), ymax = c())
  intron_line <- data.frame(x = c(), y = c(), xend = c(), yend = c())
  exon_intron_junction <- c()
  
  starts <- as.numeric(strsplit(target_gene_info$V10, split=",")[[1]])
  ends <- as.numeric(strsplit(target_gene_info$V11, split=",")[[1]])
    
  for (j in 1:length(starts)) {
      
    label <- paste(ifelse(target_gene_info$V4 == "+", j, length(starts) - j + 1))
    
    if (starts[j] <= end_target & ends[j] >= start_target) {
        
      # completely outside coding region
      if (ends[j] <= target_gene_info$V7 | starts[j] >= target_gene_info$V8) {
        txmin <- max(starts[j], start_target)
        txmax <- min(ends[j], end_target)
        tymin <- 0.4 + 1
        tymax <- 0.6 + 1
        exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
        exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 0.5, label = label))
      } else if (starts[j] >= target_gene_info$V7 & ends[j] <= target_gene_info$V8) { # completely inside coding region
        txmin <- max(starts[j], start_target)
        txmax <- min(ends[j], end_target)
        tymin <- 0.2 + 1
        tymax <- 0.8 + 1
        exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
        exon_num <- rbind(exon_num, data.frame(x = 0.5 * (txmin + txmax), y = 0.5, label = label))
      } else {
        if (starts[j] <= target_gene_info$V7) {
          txmin <- max(starts[j], start_target)
          txmax <- min(target_gene_info$V7, end_target)
          tymin <- 0.4 + 1
          tymax <- 0.6 + 1
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
            
          txmin <- max(target_gene_info$V7, start_target)
          txmax <- min(ends[j], end_target)
          tymin <- 0.2 + 1
          tymax <- 0.8 + 1
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
          exon_num <- rbind(exon_num, data.frame(
            x = 0.5 * (max(starts[j], start_target) + min(ends[j], end_target)), 
            y = 0.5, 
            label = label))
          
        } else if (ends[j] >= target_gene_info$V8) {
          txmin <- max(starts[j], start_target)
          txmax <- min(target_gene_info$V8, end_target)
          tymin <- 0.2 + 1
          tymax <- 0.8 + 1
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))
            
          txmin <- max(target_gene_info$V8, start_target)
          txmax <- min(ends[j], end_target)
          tymin <- 0.4 + 1
          tymax <- 0.6 + 1
          exon_box <- rbind(exon_box, data.frame(xmin = txmin, xmax = txmax, ymin = tymin, ymax = tymax))         
          exon_num <- rbind(exon_num, data.frame(
            x = 0.5 * (max(starts[j], start_target) + min(ends[j], end_target)), 
            y = 0.5, 
            label = label)) 
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
        ty <- 0.5 + 1
        tyend <- 0.5 + 1
        intron_line <- rbind(intron_line, data.frame(x = tx, xend = txend, y = ty, yend = tyend)) 
      }
    }
      
  }
    
  
  ##########
  # get mut info
  mut_data <- strsplit(mutation_key, ',')[[1]]
  junc_pos <- get_junction_pos(unique(sp_mut_filt$Motif_Pos), unique(sp_mut_filt$Mutation_Type))
  mut_line <- data.frame(x = junc_pos, y = 0.5 + 1)
  ##########
  
  print_info <- ggplot()
  
  if (nrow(exon_num) > 0) {
    print_info <- print_info + geom_text(data = exon_num, aes(x = x, y = y, label = label), size = 2.5)
  }

  print_info <- print_info + geom_text(data = sr_num, aes(x = x, y = y, label = label), size = 1.8)
 


  if (nrow(exon_box) > 0) {
    print_info <- print_info + geom_rect(data = exon_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey60")
  }
  
  if (nrow(intron_line) > 0) {
    print_info <- print_info + geom_segment(data = intron_line, aes(x = x, xend = xend, y = y, yend = yend), size = 0.3, color = "gray60", arrow = arrow(length = unit(0.1, "cm"))) 
  }

  example_title <- substitute(paste(a, ", ", italic(b), " (", c, ")", sep = ""), list(a = sample_name, b = gene_symbol, c = target_gene_info$V2))

  print_info <- print_info + theme_bw() +
    # ggtitle(paste(gene_symbol, " (", target_gene_info$V2, ")", sep = "")) +
    ggtitle(example_title) +
    theme(title = element_text(size = 7),
          # legend.margin = margin(0.5, 2.5, 0.5, 2.5),
          # legend.key.size = unit(1.0, "lines"), 
          # legend.key.width = unit(1.5, "lines"),
          # legend.spacing.x = unit(2, "lines"),
          legend.text = element_text(size = 6),
          legend.position = "bottom",
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_line(color="grey50", linetype="dotted"),
          panel.border=element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    # scale_x_continuous(limits = c(start_target, end_target), minor_breaks = exon_intron_junction) +
    scale_y_continuous(trans = "reverse") + 
    guides(fill=FALSE) +
    labs(color = "", linetype = "")
    # ggtitle(paste(gene_symbol, mutation_key))

  if (target_gene_info$V4 == "+") {
    print_info <- print_info + scale_x_continuous(limits = c(start_target, end_target), minor_breaks = exon_intron_junction)
  } else {
    print_info <- print_info + scale_x_continuous(trans = "reverse", limits = c(end_target, start_target), minor_breaks = exon_intron_junction)
  }
  
  print_info <- print_info + geom_point(data = mut_line, aes(x = x, y = y), color = "red", shape = 4)

  if (nrow(major_simbol) > 0) {
    print_info <- print_info + geom_point(data = major_simbol, aes(x = x, y = y), color = "#ff7f00", shape = 19)
  }

  print_info <- print_info + geom_segment(data = splicing_line, aes(x = x, xend = xend, y = y, yend = yend, 
                                                                    color = splicing_class, linetype = is_inframe), size = 0.4) +
    # scale_color_manual(values = splicing_class_colour) +
    scale_color_manual(values = splicing_class_colour[c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS", "Intron retention")]) + 
    scale_linetype_manual(values = c("In-frame" = "solid", "Frameshift" = "dashed"))

           
  return(print_info)                        
  
}



get_most_freq_splicing_calc_expression <- function(splicing_mutation, sample_name, gene_symbol, mutation_key, major_sp) {
 
  sp_mut_filt <- splicing_mutation %>% filter(Mutation_Key == mutation_key, Sample_Name == sample_name)

  # get gene info
  tmp_ref_gene <- refGene %>% filter(V13 == gene_symbol)
  gene_ind <- order(as.numeric(sub("NM_", "", tmp_ref_gene$V2)))[1]
  target_gene_info <- tmp_ref_gene[gene_ind,]

  intron_ind <- grepl("Intron retention", sp_mut_filt$Splicing_Class)

  sp_mut_filt$Splicing_Key[intron_ind] <-
    unlist(lapply(sp_mut_filt$Splicing_Key[intron_ind],
                  function(x) {get_intron_key(x, target_gene_info$V2)}))

  sp_mut_filt <- sp_mut_filt %>%
    select(Cancer_Type, Gene_Symbol, Sample_Name, Mutation_Key,
           Motif_Pos, Mutation_Type, Splicing_Key, Splicing_Class, Is_Inframe, Supporting_Read_Num) %>%
    group_by(Cancer_Type, Gene_Symbol, Sample_Name, Mutation_Key,
             Motif_Pos, Mutation_Type, Splicing_Key, Splicing_Class, Is_Inframe) %>%
    summarize(mSupporting_Read_Num = max(Supporting_Read_Num)) %>%
    arrange(desc(mSupporting_Read_Num))


  most_num <- 0
  read_num_vec <- c()
  for (i in 1:nrow(sp_mut_filt)) {
    if (as.character(sp_mut_filt[i, "Splicing_Key"]) == major_sp) most_num <- as.numeric(sp_mut_filt[i, "mSupporting_Read_Num"])
    read_num_vec <- c(read_num_vec, as.numeric(sp_mut_filt[i, "mSupporting_Read_Num"]))
  }
   
  substitute(paste(frac(a, b), " = ", c), list(a = most_num , b = paste(read_num_vec, collapse = " + "), c = round(most_num / sum(read_num_vec), digits = 4)))

}



get_legend_info <- function() {

    dummy_seg <- data.frame(x = 1:4, xend = 2:5, y = 1:4, yend = 2:5, 
                            splicing_class = c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS", "Intron retention"),
                            is_inframe = c("In-frame", "Frameshift", "In-frame", "Frameshift"))

    dummy_seg$splicing_class <- factor(dummy_seg$splicing_class,
                                       levels = c("Exon skipping", "Alternative 5'SS", "Alternative 3'SS", "Intron retention"))

    dummy_seg$is_inframe <- factor(dummy_seg$is_inframe,
                                   levels = c("In-frame", "Frameshift"))

    p <- ggplot() + 
           geom_segment(data = dummy_seg, aes(x = x, xend = xend, y = y, yend = yend, color = splicing_class, linetype = is_inframe)) +
           theme(legend.position = "bottom",
                 legend.text = element_text(size = 6)) +
           labs(color = "", linetype = "") +
           scale_color_manual(values = splicing_class_colour) +
           scale_linetype_manual(values = c("In-frame" = "solid", "Frameshift" = "dashed")) +
           guides(color = guide_legend(order = 1), linetype = guide_legend(order = 0))

    return(p)

}



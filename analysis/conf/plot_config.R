
my_theme <- function() {
   theme_bw(base_family = "Helvetica") %+replace% 
   theme(title = element_text(size = 7),
         panel.border = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "grey20", size = 0.5), 
         axis.text = element_text(size = 7),
         axis.title = element_text(size = 7),
         legend.key = element_blank(), 
         legend.key.size = unit(0.25, "cm"),
         legend.margin = margin(0.5, 0.5, 0.5, 0.5),
         legend.text = element_text(size = 6),
         strip.text = element_text(size = 6),
         # strip.background = element_rect(fill = "white", colour = "black", size = 1),
         strip.background = element_blank(), 
         complete = TRUE)
} 

splicing_class_colour <- c(
  "Exon skip" = "#fb8072",
  "Exon skipping" = "#fb8072",  
  "Alternative 5' splice site" = "#80b1d3",
  "Alternative 5'-ss" = "#80b1d3",
  "Alternative 5'SS" = "#80b1d3",
  "Alternative 3' splice site" = "#b3de69",
  "Alternative 3'-ss" = "#b3de69",
  "Alternative 3'SS" = "#b3de69",    
  "Intron retention" = "#bebada",
  "Complex" = "#fdb462",
  "Exon skip (frameshift)" = "#fb8072",
  "Exon skip (inframe)" = "#fb8072",
  "Exon skipping (frameshift)" = "#fb8072",
  "Exon skipping (inframe)" = "#fb8072",
  "Alternative 5' splice site (frameshift)" = "#80b1d3",
  "Alternative 5' splice site (inframe)" = "#80b1d3", 
  "Alternative 3' splice site (frameshift)" = "#b3de69", 
  "Alternative 3' splice site (inframe)" = "#b3de69",  
  "Alternative 5'-ss (frameshift)" = "#80b1d3",
  "Alternative 5'-ss (inframe)" = "#80b1d3", 
  "Alternative 3'-ss (frameshift)" = "#b3de69",
  "Alternative 3'-ss (inframe)" = "#b3de69",  
  "Alternative 5'SS (frameshift)" = "#80b1d3",
  "Alternative 5'SS (inframe)" = "#80b1d3",
  "Alternative 3'SS (frameshift)" = "#b3de69",
  "Alternative 3'SS (inframe)" = "#b3de69",
  "No change" = "#d9d9d9",
  "Silent" = "#d9d9d9",
  "Missense" = "#d9d9d9", 
  "Nonsense" = "#d9d9d9",
  "Inframe indel" = "#d9d9d9", 
  "Frameshift indel" = "#d9d9d9"
)

signature_colour <- c(
  "Age" = "#e41a1c", 
  "APOBEC" = "#377eb8", 
  "Tobacco" = "#4daf4a",
  "MMR defect" = "#984ea3",
  "Ultraviolet" =  "#ff7f00",
  "POLE" = "#ffff33",
  "Microsatellite" = "#a65628",
  "Other" = "#f781bf"
)

exon_size_d <- 3
intron_size_d <- 6
intron_size_a <- 6
exon_size_a <- 1


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}



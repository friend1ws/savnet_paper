library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

splicing_mut_info <- read.table("../output/omega.splicing_mutation.info.txt", header = TRUE, sep = "\t", 
                                as.is=TRUE, quote="", stringsAsFactors = FALSE)


mut_count <- read.table("../../matome/omega.mut_count.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

##########
# removing samples with extreme number of mutations ??
type2thres <- mut_count %>% group_by(Cancer_Type) %>% summarize(q90_count = quantile(Mut_Count, probs=1))

ok_list_sample <- c()
for (ctype in type2thres$Cancer_Type) {
  temp_list <- mut_count %>% 
    filter(Cancer_Type == ctype, Mut_Count <= type2thres$q90_count[type2thres$Cancer_Type == ctype]) %>% 
    select(Sample_Name)
  
  ok_list_sample <- c(ok_list_sample, temp_list$Sample_Name)
}

# splicing_mut_info_filt <- splicing_mut_info %>% filter(Sample_Name %in% ok_list_sample)
splicing_mut_info_filt <- splicing_mut_info 
##########

splicing_mut_info_filt_snv <- splicing_mut_info_filt %>% 
  filter(Ref_Mut != "-", Alt_Mut != "-")
  
##########
# organize the splicing type information
splicing_mut_info_filt_snv[
  grep(";", splicing_mut_info_filt_snv$GenomonSplicingMutation),
  "GenomonSplicingMutation"] <- "complex"

splicing_mut_info_filt_snv[
  splicing_mut_info_filt_snv$GenomonSplicingMutation == "---",
  "GenomonSplicingMutation"] <- "no-change"

splicing_mut_info_filt_snv[
  splicing_mut_info_filt_snv$GenomonSplicingMutation == "opposite-side-intron-retention",
  "GenomonSplicingMutation"] <- "intron-retention"

splicing_mut_info_filt_snv[
  splicing_mut_info_filt_snv$GenomonSplicingMutation == "intronic-alternative-5'-splice-site",
  "GenomonSplicingMutation"] <- "alternative-5'-splice-site"

splicing_mut_info_filt_snv[
  splicing_mut_info_filt_snv$GenomonSplicingMutation == "intronic-alternative-3'-splice-site",
  "GenomonSplicingMutation"] <- "alternative-3'-splice-site"

splicing_mut_info_filt_snv$GenomonSplicingMutation <-
  factor(splicing_mut_info_filt_snv$GenomonSplicingMutation,
         levels = c("exon-skip", "alternative-5'-splice-site", "alternative-3'-splice-site",
                    "intron-retention", "complex"))
##########


##########
# obtain the ratio
splicing_mut_info_filt_snv_count <- splicing_mut_info_filt_snv %>% 
  group_by(Rel_Start_Motif, Type_Motif, GenomonSplicingMutation) %>% 
  summarize(count = n())

splicing_mut_info_filt_snv_count_total <- splicing_mut_info_filt_snv_count %>% 
  group_by(Rel_Start_Motif, Type_Motif) %>% 
  summarize(total_count = sum(count))

splicing_mut_info_filt_snv_ratio <- 
  left_join(splicing_mut_info_filt_snv_count, splicing_mut_info_filt_snv_count_total, by = c("Rel_Start_Motif", "Type_Motif")) %>%
  mutate(ratio = count / total_count) 


pos_colour <- rep("grey30", 8)
pos_colour[3:4] <- "red"

p_donor <- ggplot(splicing_mut_info_filt_snv_ratio %>% 
         filter(Type_Motif == "donor", GenomonSplicingMutation != "no-change"), 
       aes(x = Rel_Start_Motif, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", fill = "") +
  ggtitle("donor") +
  ylim(c(0, 0.2)) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:8, 
                   labels = c("A", "G", "G", "T", "R", "A", "G", "T")) +
  guides(fill = FALSE)
  

  
pos_colour <- rep("grey30", 9)
pos_colour[7:8] <- "red"
  
p_acceptor <- ggplot(splicing_mut_info_filt_snv_ratio %>% 
         filter(Type_Motif == "acceptor", GenomonSplicingMutation != "no-change"), 
       aes(x = Rel_Start_Motif, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", fill = "") +
  ggtitle("acceptor") +
  ylim(c(0, 0.2)) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:9, 
                   labels = c("Y", "Y", "Y", "Y", "N", "C", "A", "G", "G")) +
  guides(fill = FALSE)




g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


p_dummy_for_legend <- ggplot(splicing_mut_info_filt_snv_ratio %>%
                               filter(GenomonSplicingMutation != "no-change"), 
                             aes(x = Rel_Start_Motif, y = ratio, fill = GenomonSplicingMutation)) + 
  geom_bar(stat = "identity") +
  labs(x = "", fill = "") +
  scale_fill_brewer(palette = "Pastel1") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow=2, byrow=TRUE))


p_donor_acceptor <- plot_grid(p_donor, p_acceptor, ncol = 2, align = "h")


plot_grid(p_donor_acceptor, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.2))

ggsave("../output/splicing_mut_count.png", width = 8, height = 4)



library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)


mem_info <- read.table("../output/omega.mut_membership.result.txt", sep = "\t", header = TRUE) %>%
  filter(Is_GSM == "TRUE") 


sig_type_order <- c("Age", "APOBEC", "Tobacco", "MMR defect", "Ultraviolet",
                    "POLE", "Microsatellite", "Other")

mem_info$Sig_Type <- "Other"
mem_info$Sig_Type[mem_info$COSM_ID == "4" & mem_info$Corr >= 0.75] <- "Tobacco"
mem_info$Sig_Type[mem_info$COSM_ID == "1" & mem_info$Corr >= 0.75] <- "Age"
mem_info$Sig_Type[mem_info$COSM_ID == "6" & mem_info$Corr >= 0.75] <- "MMR defect"
mem_info$Sig_Type[mem_info$COSM_ID %in% c("2", "13") & mem_info$Corr >= 0.75] <- "APOBEC"
mem_info$Sig_Type[mem_info$COSM_ID == "21" & mem_info$Corr >= 0.75] <- "Microsatellite"
mem_info$Sig_Type[mem_info$COSM_ID %in% c("10", "11", "32") & mem_info$Corr >= 0.75] <- "POLE"
mem_info$Sig_Type[mem_info$COSM_ID == "7" & mem_info$Corr >= 0.75] <- "Ultraviolet"

mem_info$Sig_Type <- factor(mem_info$Sig_Type, levels = sig_type_order)


pos_colour <- rep("grey30", 10)
pos_colour[4:5] <- "red"

p_donor <- ggplot(mem_info %>% filter(Splice_Site == "donor"), aes(x = Splice_Pos, y = Membership_Sum, fill = Sig_Type)) + 
  geom_bar(stat = "identity") +
  ggtitle("donor") +
  scale_fill_manual(values = c("Age" = "#e41a1c", "APOBEC" = "#377eb8", "Tobacco" = "#4daf4a", "MMR defect" = "#984ea3",
                               "Ultraviolet" =  "#ff7f00", "POLE" = "#ffff33", "Microsatellite" = "#a65628", "Other" = "#f781bf")) +
  labs(x = "", y = "Signature membership", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:9, 
                   labels = c("M", "A", "G", "G", "T", "R", "A", "G", "T")) +
  guides(fill = FALSE)


pos_colour <- rep("grey30", 7)
pos_colour[5:6] <- "red"


p_acceptor <- ggplot(mem_info %>% filter(Splice_Site == "acceptor"), aes(x = Splice_Pos, y = Membership_Sum, fill = Sig_Type)) + 
  geom_bar(stat = "identity") +
  ggtitle("acceptor") +
  scale_fill_manual(values = c("Age" = "#e41a1c", "APOBEC" = "#377eb8", "Tobacco" = "#4daf4a", "MMR defect" = "#984ea3",
                               "Ultraviolet" =  "#ff7f00", "POLE" = "#ffff33", "Microsatellite" = "#a65628", "Other" = "#f781bf")) +
  labs(x = "", y = "Signature membership", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = pos_colour, size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:7, 
                   labels = c("Y", "Y", "N", "C", "A", "G", "G")) +
  guides(fill = FALSE)



g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


p_dummy_for_legend <- ggplot(mem_info, aes(x = Splice_Pos, y = Membership_Sum, fill = Sig_Type)) + 
  geom_bar(stat = "identity") +
  labs(x = "", fill = "") +
  scale_fill_manual(values = c("Age" = "#e41a1c", "APOBEC" = "#377eb8", "Tobacco" = "#4daf4a", "MMR defect" = "#984ea3",
                               "Ultraviolet" =  "#ff7f00", "POLE" = "#ffff33", "Microsatellite" = "#a65628", "Other" = "#f781bf")) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow=2, byrow=TRUE))



p_donor_acceptor <- plot_grid(p_donor, p_acceptor, ncol = 2, rel_widths = c(1, 0.9), align = "h")

plot_grid(p_donor_acceptor, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.2))


ggsave("../output/signature_membership.pdf", width = 8, height = 4)

##########
###

mem_info_all <- read.table("../output/omega.mut_membership.all.result.txt", sep = "\t", header = TRUE)

sig_type_order <- c("Age", "APOBEC", "Tobacco", "MMR defect", "Ultraviolet",
                    "POLE", "Microsatellite", "Other")

mem_info_all$Sig_Type <- "Other"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "4" & mem_info_all$Corr >= 0.75] <- "Tobacco"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "1" & mem_info_all$Corr >= 0.75] <- "Age"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "6" & mem_info_all$Corr >= 0.75] <- "MMR defect"
mem_info_all$Sig_Type[mem_info_all$COSM_ID %in% c("2", "13") & mem_info_all$Corr >= 0.75] <- "APOBEC"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "21" & mem_info_all$Corr >= 0.75] <- "Microsatellite"
mem_info_all$Sig_Type[mem_info_all$COSM_ID %in% c("10", "31", "32") & mem_info_all$Corr >= 0.75] <- "POLE"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "7" & mem_info_all$Corr >= 0.75] <- "Ultraviolet"

mem_info_all$Sig_Type <- factor(mem_info_all$Sig_Type, levels = sig_type_order)


base_ratio <- sum(mem_info$Membership_Sum) / sum(mem_info_all$Membership_Sum)


sig2mem <- mem_info %>% 
  group_by(Sig_Type) %>% 
  summarize(mem_sum = sum(Membership_Sum)) 

sig2mem_all <- mem_info_all %>% 
  group_by(Sig_Type) %>% 
  summarize(mem_sum_all = sum(Membership_Sum)) 


gsm_ratio <- sig2mem %>% left_join(sig2mem_all, by = c("Sig_Type")) %>% 
  mutate(gsm_ratio = (mem_sum / mem_sum_all) / base_ratio)


ggplot(gsm_ratio %>% filter(Sig_Type != "Other"), aes(x = Sig_Type, y = gsm_ratio, fill = Sig_Type)) + 
  geom_bar(stat = "identity") + coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Relative splicing mutation ratio", fill = "") +
  theme(legend.position = "bottom") +
  scale_x_discrete(limits = rev(sig_type_order[1:7])) +
  scale_fill_manual(values = c("Age" = "#e41a1c", "APOBEC" = "#377eb8", "Tobacco" = "#4daf4a", "MMR defect" = "#984ea3",
                               "Ultraviolet" =  "#ff7f00", "POLE" = "#ffff33", "Microsatellite" = "#a65628", "Other" = "#f781bf"))

 
ggsave("../output/rel_sp_ratio.pdf", width = 6, height = 3)



##########
# LUSC

p_ctype <-  list(NA, NA, NA, NA, NA, NA)
type_num <- rep(0, 6)
ctype_vec <- c("LUSC", "LUAD", "HNSC", "SKCM", "COAD", "UCEC")

for(i in 1:length(ctype_vec)) {
  
  ctype <- ctype_vec[i]

  sig2mem <- mem_info %>% filter(Cancer_Type == ctype) %>%
    group_by(Sig_Type) %>% 
    summarize(mem_sum = sum(Membership_Sum)) 

  sig2mem_all <- mem_info_all %>% filter(Cancer_Type == ctype) %>%
    group_by(Sig_Type) %>% 
    summarize(mem_sum_all = sum(Membership_Sum)) 


  gsm_ratio <- sig2mem %>% left_join(sig2mem_all, by = c("Sig_Type")) %>% 
    mutate(gsm_ratio = (mem_sum / mem_sum_all) / base_ratio)



  tp <- ggplot(gsm_ratio, aes(x = Sig_Type, y = gsm_ratio, fill = Sig_Type)) + 
    geom_bar(stat = "identity") + coord_flip() +
    theme_minimal() +
    ggtitle(ctype) + 
    # scale_fill_brewer(palette = "Set1") +
    labs(x = "", y = "Relative splicing mutation ratio", fill = "") +
    theme(legend.position = "bottom") +
    scale_x_discrete(limits = rev(sig_type_order[sig_type_order %in% gsm_ratio$Sig_Type])) +
    scale_fill_manual(values = c("Age" = "#e41a1c", "APOBEC" = "#377eb8", "Tobacco" = "#4daf4a", "MMR defect" = "#984ea3",
                                 "Ultraviolet" =  "#ff7f00", "POLE" = "#ffff33", "Microsatellite" = "#a65628", "Other" = "#f781bf")) +
    guides(fill = FALSE)

  p_ctype[[i]] <- tp
  type_num[i] <- nrow(gsm_ratio)
}



p_tobacco <- plot_grid(p_ctype[[1]], p_ctype[[2]], nrow = 2)
p_uv <- plot_grid(p_ctype[[3]], p_ctype[[4]], nrow = 2)
p_pole <- plot_grid(p_ctype[[5]], p_ctype[[6]], nrow = 2)

plot_grid(plot_grid(p_tobacco, p_uv, p_pole, nrow = 1), g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.1))

ggsave(paste("../output/rel_sp_ratio_ctype.pdf", sep = ""), width = 10, height = 7.5)





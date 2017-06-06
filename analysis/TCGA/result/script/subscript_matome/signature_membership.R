library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(scales)

source("../../../conf/plot_config.R")

mem_info <- read.table("../temporary/TCGA.mut_membership.savnet.result.txt", sep = "\t", header = TRUE) %>%
  filter(Is_GSM == "TRUE") 


sig_type_order <- c("Aging", "APOBEC", "Smoking", "MMR defect", "Ultraviolet",
                    "POLE", "Microsatellite", "Other")


mem_info$Sig_Type <- "Other"
mem_info$Sig_Type[mem_info$COSM_ID == "4" & mem_info$Corr >= 0.75] <- "Smoking"
mem_info$Sig_Type[mem_info$COSM_ID == "1" & mem_info$Corr >= 0.75] <- "Aging"
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
  ggtitle("Donor") +
  scale_fill_manual(values = signature_colour) +
  labs(x = "", y = "Signature membership", fill = "") +
  my_theme() +
  theme(axis.text.x = element_text(colour = pos_colour), # , size = rel(1.5)),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:9, 
                   labels = c("M", "A", "G", "G", "T", "R", "A", "G", "T")) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = FALSE)


pos_colour <- rep("grey30", 7)
pos_colour[5:6] <- "red"


p_acceptor <- ggplot(mem_info %>% filter(Splice_Site == "acceptor"), aes(x = Splice_Pos, y = Membership_Sum, fill = Sig_Type)) + 
  geom_bar(stat = "identity") +
  ggtitle("Acceptor") +
  scale_fill_manual(values = signature_colour) +
  labs(x = "", y = "Signature membership", fill = "") +
  my_theme() +
  theme(axis.text.x = element_text(colour = pos_colour), # , size = rel(1.5)),
        # axis.text.y = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.2)),
        # legend.text = element_text(size = rel(1)),
        # legend.title = element_text(size = rel(1)),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(limits = 1:7, 
                   labels = c("Y", "Y", "N", "C", "A", "G", "G")) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = FALSE)



g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


p_dummy_for_legend <- ggplot(mem_info, aes(x = Splice_Pos, y = Membership_Sum, fill = Sig_Type)) + 
  geom_bar(stat = "identity") +
  labs(x = "", fill = "") +
  scale_fill_manual(values = signature_colour) +
  my_theme() +
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(nrow=2, byrow=TRUE))



p_donor_acceptor <- plot_grid(p_donor, p_acceptor, ncol = 2, rel_widths = c(1, 0.9), align = "h")

plot_grid(p_donor_acceptor, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.1))


ggsave("../figure/signature_membership.tiff", width = 15, height = 6, dpi = 600, units = "cm")

##########
###

mem_info_all <- read.table("../temporary/TCGA.mut_membership.all.result.txt", sep = "\t", header = TRUE)

sig_type_order <- c("Aging", "APOBEC", "Smoking", "MMR defect", "Ultraviolet",
                    "POLE", "Microsatellite", "Other")
sig_type_order2 <- c("Aging", "APOBEC", "Smoking", "Ultraviolet", "POLE")


mem_info_all$Sig_Type <- "Other"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "4" & mem_info_all$Corr >= 0.75] <- "Smoking"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "1" & mem_info_all$Corr >= 0.75] <- "Aging"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "6" & mem_info_all$Corr >= 0.75] <- "MMR defect"
mem_info_all$Sig_Type[mem_info_all$COSM_ID %in% c("2", "13") & mem_info_all$Corr >= 0.75] <- "APOBEC"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "21" & mem_info_all$Corr >= 0.75] <- "Microsatellite"
mem_info_all$Sig_Type[mem_info_all$COSM_ID %in% c("10", "31", "32") & mem_info_all$Corr >= 0.75] <- "POLE"
mem_info_all$Sig_Type[mem_info_all$COSM_ID == "7" & mem_info_all$Corr >= 0.75] <- "Ultraviolet"

mem_info_all$Sig_Type <- factor(mem_info_all$Sig_Type, levels = sig_type_order)


base_ratio <- sum(mem_info$Membership_Sum) / sum(mem_info_all$Membership_Sum)
# base_ratio <- 1

sig2mem <- mem_info %>% 
  group_by(Sig_Type) %>% 
  summarize(mem_sum = sum(Membership_Sum)) 

sig2mem_all <- mem_info_all %>% 
  group_by(Sig_Type) %>% 
  summarize(mem_sum_all = sum(Membership_Sum)) 


gsm_ratio <- sig2mem %>% left_join(sig2mem_all, by = c("Sig_Type")) %>% 
  mutate(gsm_ratio = (mem_sum / mem_sum_all)) 


ggplot(gsm_ratio %>% filter(Sig_Type %in% sig_type_order2), aes(x = Sig_Type, y = gsm_ratio, fill = Sig_Type)) + 
  geom_bar(stat = "identity") + # coord_flip() +
  geom_abline(intercept = base_ratio, slope = 0, colour = "#a65628", alpha = 0.9, linetype = "longdash") + 
  my_theme() +
  labs(x = "", y = "SAV fraction", fill = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_discrete(limits = sig_type_order2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = signature_colour) +
  guides(fill = FALSE)

 
ggsave("../figure/rel_sp_ratio.tiff", width = 4, height = 7, dpi = 600, units = "cm")



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

  base_ratio <- sum(sig2mem$mem_sum) / sum(sig2mem_all$mem_sum_all)


  gsm_ratio <- sig2mem %>% left_join(sig2mem_all, by = c("Sig_Type")) %>% 
    mutate(gsm_ratio = (mem_sum / mem_sum_all))



  tp <- ggplot(gsm_ratio %>% filter(Sig_Type %in% sig_type_order2), 
               aes(x = Sig_Type, y = gsm_ratio, fill = Sig_Type)) + 
    geom_bar(stat = "identity") + # coord_flip() +
    geom_abline(intercept = base_ratio, slope = 0, colour = "#a65628", alpha = 0.9, linetype = "longdash") +
    my_theme() +
    ggtitle(ctype) + 
    # scale_fill_brewer(palette = "Set1") +
    labs(x = "", y = "", fill = "") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_discrete(limits = sig_type_order2[sig_type_order2 %in% gsm_ratio$Sig_Type]) +
    scale_y_continuous(expand = c(0, 0), # labels = scientific, 
                       # breaks = seq(0, max(gsm_ratio$gsm_ratio), 0.002),
                       # minor_breaks = seq(0.001, max(gsm_ratio$gsm_ratio), 0.002),
                       limits = c(0, 0.006)) +
    scale_fill_manual(values = signature_colour) +
    guides(fill = FALSE)

  p_ctype[[i]] <- tp
  type_num[i] <- nrow(gsm_ratio)
}



# p_tobacco <- plot_grid(p_ctype[[1]], p_ctype[[2]], nrow = 2, align = "v")
# p_uv <- plot_grid(p_ctype[[3]], p_ctype[[4]], nrow = 2, align = "v")
# p_pole <- plot_grid(p_ctype[[5]], p_ctype[[6]], nrow = 2, align = "v")

#  plot_grid(p_tobacco, p_uv, p_pole, nrow = 1, align = "h")

ytitle <- ggdraw() + draw_label("SAV fraction", size = 8, angle = 90)

plot_grid(ytitle, plot_grid(plotlist = p_ctype, nrow = 2, align = "hv"), nrow = 1, rel_widths = c(0.1, 1), scale = 0.99)

ggsave(paste("../figure/rel_sp_ratio_ctype.tiff", sep = ""), width = 10, height = 9, dpi = 600, units = "cm")


##########
library(pmsignature)


my_sig_num <- read.table("../../data/pmsignature/db/manual_sig_num.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE)

get_plot <- function(sig_num, title) {
  selected_key <- mem_info_all %>% filter(COSM_ID == as.character(sig_num)) %>% filter(Corr == max(Corr)) 
  sig_num <- my_sig_num[my_sig_num[,1] == selected_key$Cancer_Type, 2]
  load(paste("../../data/pmsignature/", selected_key$Cancer_Type, "/pmsignature/ind.", sig_num, ".Rdata", sep = ""))
  visPMSignature(resultForSave[[1]], selected_key$Sig_Num, charSize = 2) + 
	ggtitle(title) + 
	theme(panel.background = element_blank(),
	      title = element_text(size = 7))
}

p_age <- get_plot(1, "Aging")
p_appobec <- get_plot(13, "APOBEC")
p_tobacco <- get_plot(4, "Smoking")
# p_mmr <- get_plot(6, "MMR defect")
p_uv <- get_plot(7, "Ultraviolet")
p_pole1 <- get_plot(31, "POLE1")
p_pole2 <- get_plot(32, "POLE2")
# p_ms <- get_plot(21, "Microsatellite")


theme_set(theme_gray())

# plot_grid(p_age, p_appobec, p_tobacco, p_mmr, p_uv, p_pole1, p_pole2, p_ms, nrow = 2)
plot_grid(p_age, p_appobec, p_tobacco, p_uv, p_pole1, p_pole2, nrow = 3)

ggsave("../figure/pmsignature_list.tiff", width = 8, height = 9, dpi = 600, units = "cm") 



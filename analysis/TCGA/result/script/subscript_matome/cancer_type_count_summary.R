library(dplyr)
library(ggplot2)
library(ggrepel)

source("../../../conf/plot_config.R")

mut_count <- read.table("../temporary/omega.mut_count.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


sv_count_median <- mut_count %>% 
  filter(SAV_Count <= 8) %>% 
  group_by(Cancer_Type) %>% 
  summarize(median_count = median(SAV_Count), q25_count = quantile(SAV_Count, probs=0.25), q75_count = quantile(SAV_Count, probs=0.75)) %>% 
  arrange(median_count, q75_count, q25_count, Cancer_Type)

mut_count$Cancer_Type <- factor(mut_count$Cancer_Type, levels = sv_count_median$Cancer_Type)


ggplot(mut_count, aes(x = Cancer_Type, y = SAV_Count)) + 
  geom_point(position = position_jitter(width = 0.1, height = 0.1), colour = "#fc8d59", alpha = 0.2, size = 0.3) +
  geom_boxplot(outlier.colour = NA, fill = "#3288bd", size = 0.3) +
  my_theme() +
  labs(y = "SAV count", x = "Cancer type") +
  theme(# axis.text = element_text(size = rel(1.2)),
        # axis.title = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()) + # ,
        # legend.text = element_text(size = rel(1.2)),
        # legend.title = element_text(size = rel(1.5))) +
  scale_y_continuous(limit = c(-0.2, 8))

ggsave("../figure/count_summary.tiff", width = 9, height = 4, dpi = 600, units = "cm")

##########

mut_count_median <- mut_count %>%
  mutate(Mut_Count = SNV_Count + Indel_Count) %>%
  group_by(Cancer_Type) %>%
  summarize(SAV_Count_Median = median(SAV_Count), Mut_Count_Median = median(Mut_Count))

ggplot(mut_count_median) +
  geom_point(aes(x = Mut_Count_Median, y = SAV_Count_Median), colour = "#1a9850", size = 0.5) + 
  geom_text_repel(data = mut_count_median %>% filter(SAV_Count_Median > 0), 
                  aes(x = Mut_Count_Median, y = SAV_Count_Median, label = Cancer_Type), size = 1.8) +
  labs(x = "Median mutation count", y = "Median SAV count") +
  my_theme() # +
  # theme(axis.text = element_text(size = rel(1.2)),
  #       axis.title = element_text(size = rel(1.5))) # +
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0))

ggsave("../figure/mut2spmut_median.tiff", width = 6, height = 6, dpi = 600, units = "cm")

##########

mut_count_trmean <- mut_count %>%
  mutate(Mut_Count = SNV_Count + Indel_Count) %>%
  group_by(Cancer_Type) %>%
  summarize(SAV_Count_Mean = mean(SAV_Count, trim = 0.1), Mut_Count_Mean = mean(Mut_Count, trim = 0.1))


lm_res <- lm(SAV_Count_Mean ~ Mut_Count_Mean, 
             mut_count_trmean %>% 
            filter(!Cancer_Type %in% c("COAD", "SKCM", "UCEC")))


ggplot(mut_count_trmean) +
  geom_point(aes(x = Mut_Count_Mean, y = SAV_Count_Mean), colour = "#1a9850", size = 0.5) + 
  geom_abline(intercept = lm_res$coefficients[1], slope = lm_res$coefficients[2],
              colour = "#d73027", alpha = 0.6) +
  geom_text_repel(data = mut_count_trmean %>% filter(SAV_Count_Mean > 0.5), 
                  aes(x = Mut_Count_Mean, y = SAV_Count_Mean, label = Cancer_Type), size = 1.8) +
  labs(x = "Trucated mean of variant count", y = "Truncated mean of SAV count") +
  my_theme() +
  scale_x_continuous(breaks=seq(0, 1250, 250),limits=c(0,1300))
  # theme(axis.text = element_text(size = rel(1.2)),
  #       axis.title = element_text(size = rel(1.5))) # +
  # scale_x_continuous(expand = c(0, 0)) +
  # scale_y_continuous(expand = c(0, 0))

ggsave("../figure/mut2spmut_trmean.tiff", width = 6, height = 6, dpi = 600, units = "cm")


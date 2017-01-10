library(dplyr)
library(ggplot2)
library(ggrepel)

splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)
mut_count <- read.table("../matome/omega.mut_count.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

sp_mut_count <- splicing_mutation %>% 
  select(Cancer_Type, Sample_Name, Mutation_Key) %>% 
  distinct() %>% 
  group_by(Cancer_Type, Sample_Name) %>% 
  summarize(Mut_Sp_Count = n())

count_summary <- full_join(sp_mut_count, mut_count, by = c("Cancer_Type", "Sample_Name"))
count_summary$Mut_Sp_Count[is.na(count_summary$Mut_Sp_Count)] <- 0

count_summary <- count_summary %>% filter(!(Cancer_Type %in% c("LAML", "OV")))

sv_count_median <- count_summary %>% 
  group_by(Cancer_Type) %>% 
  summarize(median_count = median(Mut_Sp_Count), q25_count = quantile(Mut_Sp_Count, probs=0.25), q75_count = quantile(Mut_Sp_Count, probs=0.75)) %>% 
  arrange(median_count, q75_count, Cancer_Type)

count_summary$Cancer_Type <- factor(count_summary$Cancer_Type, levels = sv_count_median$Cancer_Type)


ggplot(count_summary, aes(x = Cancer_Type, y = Mut_Sp_Count)) + 
  geom_point(position = position_jitter(width = 0.4, height = 0.3), colour = "#fc8d59", alpha = 0.2) +
  geom_boxplot(outlier.colour = NA, fill = "#3288bd") +
  theme_minimal() +
  ylim(c(-0.2, 8)) +
  labs(y = "#SASMs", x = "Cancer type") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))

ggsave("../matome/count_summary.pdf", width = 10, height = 6)

##########

count_summary_median <- count_summary %>% 
  group_by(Cancer_Type) %>%
  summarize(Mut_Sp_Count_Median = median(Mut_Sp_Count), Mut_Count_Median = median(Mut_Count))

ggplot(count_summary_median) +
  geom_point(aes(x = Mut_Count_Median, y = Mut_Sp_Count_Median), colour = "red") + 
  geom_text_repel(data = count_summary_median %>% filter(Mut_Sp_Count_Median > 0), 
                  aes(x = Mut_Count_Median, y = Mut_Sp_Count_Median, label = Cancer_Type)) +
  labs(x = "Median #mutation", y = "Median #splicing_mutation") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)))

ggsave("../matome/mut2spmut_median.png", width = 8, height = 8)

##########

count_summary_trmean <- count_summary %>% 
  group_by(Cancer_Type) %>%
  summarize(Mut_Sp_Count_Mean = mean(Mut_Sp_Count, trim = 0.1), Mut_Count_Mean = mean(Mut_Count, trim = 0.1))

ggplot(count_summary_trmean) +
  geom_point(aes(x = Mut_Count_Mean, y = Mut_Sp_Count_Mean), colour = "red") + 
  geom_text_repel(data = count_summary_trmean, 
                  aes(x = Mut_Count_Mean, y = Mut_Sp_Count_Mean, label = Cancer_Type)) +
  labs(x = "Trucated mean #mutation", y = "Truncated mean #splicing_mutation") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)))

ggsave("../matome/mut2spmut_trmean.png", width = 8, height = 8)


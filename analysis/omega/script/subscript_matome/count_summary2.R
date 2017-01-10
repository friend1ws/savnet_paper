library(dplyr)
library(ggplot2)
library(ggrepel)

# splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)
mut_count <- read.table("../matome/omega.mut_count.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# sp_mut_count <- splicing_mutation %>% 
#   select(Cancer_Type, Sample_Name, Mutation_Key) %>% 
#   distinct() %>% 
#   group_by(Cancer_Type, Sample_Name) %>% 
#   summarize(Mut_Sp_Count = n())

# count_summary <- full_join(sp_mut_count, mut_count, by = c("Cancer_Type", "Sample_Name"))
# count_summary$Mut_Sp_Count[is.na(count_summary$Mut_Sp_Count)] <- 0

# count_summary <- count_summary %>% filter(!(Cancer_Type %in% c("LAML", "OV")))

sv_count_median <- mut_count %>% 
  filter(SASM_Count <= 8) %>% 
  group_by(Cancer_Type) %>% 
  summarize(median_count = median(SASM_Count), q25_count = quantile(SASM_Count, probs=0.25), q75_count = quantile(SASM_Count, probs=0.75)) %>% 
  arrange(median_count, q75_count, q25_count, Cancer_Type)

mut_count$Cancer_Type <- factor(mut_count$Cancer_Type, levels = sv_count_median$Cancer_Type)


ggplot(mut_count, aes(x = Cancer_Type, y = SASM_Count)) + 
  geom_point(position = position_jitter(width = 0.4, height = 0.3), colour = "#fc8d59", alpha = 0.2) +
  geom_boxplot(outlier.colour = NA, fill = "#3288bd") +
  theme_minimal() +
  ylim(c(-0.2, 8)) +
  labs(y = "SASM Count", x = "Cancer type") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))

ggsave("../matome/count_summary.pdf", width = 10, height = 6)

##########

mut_count_median <- mut_count %>%
  mutate(Mut_Count = SNV_Count + Indel_Count) %>%
  group_by(Cancer_Type) %>%
  summarize(SASM_Count_Median = median(SASM_Count), Mut_Count_Median = median(Mut_Count))

ggplot(mut_count_median) +
  geom_point(aes(x = Mut_Count_Median, y = SASM_Count_Median), colour = "red") + 
  geom_text_repel(data = mut_count_median %>% filter(SASM_Count_Median > 0), 
                  aes(x = Mut_Count_Median, y = SASM_Count_Median, label = Cancer_Type)) +
  labs(x = "Median mutation count", y = "Median SASM count") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)))

ggsave("../matome/mut2spmut_median.pdf", width = 8, height = 8)

##########

mut_count_trmean <- mut_count %>%
  mutate(Mut_Count = SNV_Count + Indel_Count) %>%
  group_by(Cancer_Type) %>%
  summarize(SASM_Count_Mean = mean(SASM_Count, trim = 0.1), Mut_Count_Mean = mean(Mut_Count, trim = 0.1))


ggplot(mut_count_trmean) +
  geom_point(aes(x = Mut_Count_Mean, y = SASM_Count_Mean), colour = "red") + 
  geom_text_repel(data = mut_count_trmean, 
                  aes(x = Mut_Count_Mean, y = SASM_Count_Mean, label = Cancer_Type)) +
  labs(x = "Trucated mean mutation count", y = "Truncated mean SASM count") +
  theme_bw() +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)))

ggsave("../matome/mut2spmut_trmean.pdf", width = 8, height = 8)


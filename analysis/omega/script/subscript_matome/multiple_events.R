
library(dplyr)
library(ggplot2)

splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

multiple_effect <- splicing_mutation %>% 
  group_by(Cancer_Type, Sample_Name, Gene_Symbol, Mutation_Key) %>% 
  summarize(splice_count = n()) %>%
  arrange(desc(splice_count))

multiple_effect_count <- multiple_effect %>%
  group_by(splice_count) %>% 
  summarize(count = n()) 


ggplot(multiple_effect_count, aes(x = splice_count, y = count)) + 
  geom_bar(stat = "identity") + 
  labs(x = "#splicing event", y = "#mutation") +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)))

ggsave("../matome/multi_event.png", width = 6, height = 6)


splicing_multiple_mut <- splicing_mutation %>% 
  select(Mutation_Key, Splicing_Key, Gene_Symbol, Splicing_Class, Is_Inframe) %>% 
  distinct() %>%
  group_by(Splicing_Key, Gene_Symbol, Splicing_Class, Is_Inframe) %>%
  summarize(mut_count = n()) 

splicing_multiple_mut_count <- splicing_multiple_mut %>%
  group_by(mut_count) %>% 
  summarize(count = n())


ggplot(splicing_multiple_mut_count, aes(x = mut_count, y = count)) + 
  geom_bar(stat = "identity") + 
  labs(x = "#associated mutations", y = "#splicing events") +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)))

ggsave("../matome/multi_event_sp.png", width = 6, height = 6)


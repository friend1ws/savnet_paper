library(dplyr)
library(ggplot2)

source("../../../conf/plot_config.R")

splicing_mutation <- read.table("../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt", 
                                header = TRUE, sep = "\t", quote="", stringsAsFactors = FALSE) %>%
    filter(IR_filtered == "FALSE")

multiple_effect <- splicing_mutation %>% 
  group_by(Cancer_Type, Sample_Name, Gene_Symbol, Mutation_Key) %>% 
  summarize(splice_count = n()) %>%
  arrange(desc(splice_count))

multiple_effect_count <- multiple_effect %>%
  group_by(splice_count) %>% 
  summarize(count = n()) 

print(multiple_effect_count)

multiple_effect_count$splice_count2 <- factor(multiple_effect_count$splice_count, levels = 1:max(multiple_effect_count$splice_count))

ggplot(multiple_effect_count, aes(x = splice_count2, y = log10(count + 1))) + 
  geom_bar(stat = "identity", fill = "#6baed6") + 
  labs(x = "Associated splicing event count", y = "log10(mutation count + 1)") +
  # scale_y_log10() +
  my_theme() +
  scale_y_continuous(expand = c(0, 0))

ggsave("../figure/multi_event.pdf", width = 3, height = 3)


# splicing_multiple_mut <- splicing_mutation %>% 
#   select(Mutation_Key, Splicing_Key, Gene_Symbol, Splicing_Class, Is_Inframe) %>% 
#   distinct() %>%
#   group_by(Splicing_Key, Gene_Symbol, Splicing_Class, Is_Inframe) %>%
#   summarize(mut_count = n()) 

# splicing_multiple_mut_count <- splicing_multiple_mut %>%
#   group_by(mut_count) %>% 
#   summarize(count = n())


# ggplot(splicing_multiple_mut_count, aes(x = mut_count, y = count)) + 
#   geom_bar(stat = "identity") + 
#   labs(x = "Associated mutation count", y = "log10(splicing event count)") +
#   scale_y_log10() +
#   my_theme() +
#   theme(axis.text = element_text(size = rel(1.2)),
#         axis.title = element_text(size = rel(1.2)))

# ggsave("../matome/multi_event_sp.png", width = 6, height = 6)


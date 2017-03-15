library(dplyr)

splicing_mutation <- read.table("../../output/rescue/TCGA.savnet.with_rescued.result.txt", 
    header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE) %>% 
    filter(IR_filtered != "TRUE")


splicing_mutation_count <- splicing_mutation %>% 
  select(Sample_Name, Gene_Symbol) %>%
  distinct() %>%
  group_by(Gene_Symbol) %>%
  summarize(SAV_Count = n()) %>% 
  arrange(desc(SAV_Count))

write.table(splicing_mutation_count, "../table/sav.count.txt", sep = "\t", quote = FALSE, row.names = FALSE)



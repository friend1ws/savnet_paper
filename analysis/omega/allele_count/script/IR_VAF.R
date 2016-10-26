library(ggplot2)
library(dplyr)

allele_count <- read.table("../output/omega.genomon_splicing_mutation.allele_count.txt", header = TRUE, sep = "\t", 
                           as.is=TRUE, quote="", stringsAsFactors = FALSE)

allele_count <- allele_count %>% 
  mutate(ratio = Intron_Retention_Positive / (Intron_Retention_Negative + Intron_Retention_Positive))

allele_count$Mutation_Type2 <- allele_count$Mutation_Type
allele_count$Mutation_Type[allele_count$Mutation_Type == "splicing acceptor disruption"] <- "acceptor"
allele_count$Mutation_Type[allele_count$Mutation_Type == "splicing donor disruption"] <- "donor"


a <- hist(allele_count$ratio, breaks = seq(0, 1, 0.05))
write.table(data.frame(mids = a$mids, counts = a$counts), "../output/IR_VAF.txt", quote = FALSE, row.names = FALSE, sep = "\t")

ggplot(allele_count %>% filter(! is.na(ratio)),
       aes(x = ratio, fill = Mutation_Type2)) +
  geom_histogram(bins = 50, alpha = 0.5, position="stack") +
  labs(x = "Variant Allele Frequency", fill = "disrupted splicing motif") +
  ggtitle("Intron Retention Variant Allele Frequency") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_discrete(labels = c("donor", "acceptor"))


ggsave("../output/IR_VAF.png", width = 6, height = 4.5)



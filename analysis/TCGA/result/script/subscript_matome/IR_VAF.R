library(ggplot2)
library(dplyr)

source("../../../conf/plot_config.R")

allele_count <- read.table("../temporary/TCGA.savnet.allele_count.txt", header = TRUE, sep = "\t", 
                           as.is=TRUE, quote="", stringsAsFactors = FALSE)

allele_count <- allele_count %>% 
  mutate(ratio = Intron_Retention_Positive / (Intron_Retention_Negative + Intron_Retention_Positive))

allele_count$Mutation_Type2 <- allele_count$Mutation_Type
allele_count$Mutation_Type2[allele_count$Mutation_Type == "splicing acceptor disruption"] <- "Acceptor"
allele_count$Mutation_Type2[allele_count$Mutation_Type == "splicing donor disruption"] <- "Donor"
allele_count$Mutation_Type2 = factor(allele_count$Mutation_Type2, levels = c("Donor", "Acceptor"))

a <- hist(allele_count$ratio, breaks = seq(0, 1, 0.05))
write.table(data.frame(mids = a$mids, counts = a$counts), "../table/IR_VAF.txt", quote = FALSE, row.names = FALSE, sep = "\t")

ggplot(allele_count %>% filter(! is.na(ratio)),
       aes(x = ratio)) +
  geom_histogram(bins = 50, position="stack", fill = "#2166ac") +
  labs(x = "Variant Allele Frequency", y = "Frequency") +
  # ggtitle("Intron Retention Variant Allele Frequency") +
  facet_grid(.~Mutation_Type2, scales = "free_y") +
  my_theme() +
  theme(strip.text.x = element_text(angle = 0, hjust = 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = FALSE)


ggsave("../figure/IR_VAF.pdf", width = 6, height = 2.5)



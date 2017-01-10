library(dplyr)
library(ggplot2)

omega_gsm <- read.table("../matome/omega.genomon_splicing_mutation.result.txt",
sep = "\t", header = TRUE, quote = "")

omega_FDR <- omega_gsm %>% group_by(Cancer_Type) %>% summarize(FDR =
max(Q_Value))

ggplot(omega_FDR, aes(x = Cancer_Type, y = FDR)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


ggsave("../matome/cancer_type_fdr.pdf", width = 10, height = 5)


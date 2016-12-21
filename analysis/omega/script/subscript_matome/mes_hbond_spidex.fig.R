library(dplyr)
library(ggplot2)

snv_info <- read.table("../matome/omega.motif_summary.mes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

snv_info$Mutation_Type2 <- factor(snv_info$Mutation_Type,
                                  levels = rev(c("splicing donor disruption", "splicing donor creation",
                                             "splicing acceptor disruption", "splicing acceptor creation")),
                                  labels = rev(c("Donor disruption", "Donor creation",
                                             "Acceptor disruption", "Acceptor creation")))


ggplot(snv_info, aes(x = Mutation_Type2, y = mes_mt - mes_wt)) + 
  geom_boxplot() + 
  coord_flip() +
  ylim(c(-20, 20)) +
  labs(x = "", y = "Diff. of MaxEnt score") +
  theme_bw()
  
ggsave("../matome/diff_maxent.png", width = 6, height = 3)

ggplot(snv_info %>% filter(SPIDEX != "---"), aes(x = Mutation_Type2, y = as.numeric(SPIDEX))) + 
  geom_boxplot() + 
  coord_flip() +
  labs(x = "", y = "SPIDEX") +
  theme_bw()

ggsave("../matome/spidex.png", width = 6, height = 3)


ggplot(snv_info %>% filter(Mutation_Type2 %in% c("Donor disruption", "Donor creation")),
       aes(x = Mutation_Type2, y = hbond_mt - hbond_wt)) + 
  geom_boxplot() + 
  coord_flip() +
  labs(x = "", y = "Diff. of HBond scores") +
  theme_bw()

ggsave("../matome/diff_hbond.png", width = 6, height = 1.8)




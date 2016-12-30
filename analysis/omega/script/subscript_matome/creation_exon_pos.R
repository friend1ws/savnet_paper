library(dplyr)
library(ggplot2)

a <- read.table("../matome/omega.creation_exon.txt", sep = "\t", header = TRUE, quote = "") %>%
  mutate(Pos_Ratio = (Mut_Pos - Exon_Start) / (Exon_End - Exon_Start))

Pos_Ratio2 <- a$Pos_Ratio

Pos_Ratio2[a$Exon_Strand == "-"] <- 1 - Pos_Ratio2[a$Exon_Strand == "-"] 

a$Pos_Ratio2 <- Pos_Ratio2

# a <- a %>% filter(Pos_Ratio2 >= 0) %>% filter(Pos_Ratio2 <= 1)
a_donor <- a %>% filter(Mutation_Type == "splicing donor creation")
a_acceptor <- a %>% filter(Mutation_Type == "splicing acceptor creation")

ggplot(a , aes(x = Pos_Ratio2, fill = Mutation_Type)) + 
  geom_histogram(binwidth = 0.05, colour = "gray30") +
  facet_grid(.~Mutation_Type, scales = "free") +
  guides(fill = FALSE) +
  labs(x = "Exon relative position")
       
ggsave("../matome/creation_exon_pos.png", width = 8, height = 4)


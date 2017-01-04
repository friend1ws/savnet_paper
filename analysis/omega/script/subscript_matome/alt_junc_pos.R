library(dplyr)
library(ggplot2)

a <- read.table("../matome/omega.alt_junc.txt", sep = "\t", header = TRUE, quote = "") # %>%
#   filter(Mutation_Type %in% c("splicing donor creation", "splicing acceptor creation"))


# %>% mutate(Pos_Ratio = (Mut_Pos - Exon_Start) / (Exon_End - Exon_Start))

Pos_Diff <- rep(NA, nrow(a))

ind <- a$Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") & a$Exon_Strand == "+"
Pos_Diff[ind] <- a$Junc_Pos[ind] - a$Exon_Start[ind]

ind <- a$Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation") & a$Exon_Strand == "-"
Pos_Diff[ind] <- a$Exon_End[ind] - a$Junc_Pos[ind]

ind <- a$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation") & a$Exon_Strand == "+"
Pos_Diff[ind] <- a$Exon_End[ind] - a$Junc_Pos[ind]

ind <- a$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation") & a$Exon_Strand == "-"
Pos_Diff[ind] <- a$Junc_Pos[ind] - a$Exon_Start[ind]




# Pos_Diff[a$Exon_Strand == "-"] <- a$Exon_End[a$Exon_Strand == "-"] - a$Mut_Pos[a$Exon_Strand == "-"] 


# Pos_Ratio2 <- a$Pos_Ratio

# Pos_Ratio2[a$Exon_Strand == "-"] <- 1 - Pos_Ratio2[a$Exon_Strand == "-"] 

a$Pos_Diff <- Pos_Diff



a <- a %>% filter(Pos_Diff >= -100) %>% filter(Pos_Diff <= 300)
a_donor <- a %>% filter(Mutation_Type == "splicing donor creation")
a_acceptor <- a %>% filter(Mutation_Type == "splicing acceptor creation")

a$Mutation_Type <- factor(a$Mutation_Type, 
                          levels = c("splicing donor creation", "splicing acceptor creation",
                                     "splicing donor disruption", "splicing acceptor disruption"))

ggplot(a , aes(x = Pos_Diff, fill = Mutation_Type)) + 
  geom_histogram(binwidth = 5, colour = "gray30") +
  xlim(c(-100, 300)) + 
  # facet_grid(.~Mutation_Type, scales = "free") +
  facet_wrap( ~ Mutation_Type, scales = "free", ncol = 2) +
  guides(fill = FALSE) +
  labs(x = "Exonic position")

ggsave("../matome/alt_exon_pos.png", width = 6, height = 6)


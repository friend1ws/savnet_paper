library(dplyr)
library(ggplot2)

##########
# caterogy summary plot
splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

splicing_mutation$Mutation_Type <- splicing_mutation$Mutation_Type

splicing_mutation[
  grep("splicing donor disruption", splicing_mutation$Mutation_Type),
  "Mutation_Type"] <- "splicing donor disruption"

splicing_mutation[
  grep("splicing acceptor disruption", splicing_mutation$Mutation_Type),
  "Mutation_Type"] <- "splicing acceptor disruption"

splicing_mutation[
  grep("splicing donor creation", splicing_mutation$Mutation_Type),
  "Mutation_Type"] <- "splicing donor creation"

splicing_mutation[
  grep("splicing acceptor creation", splicing_mutation$Mutation_Type),
  "Mutation_Type"] <- "splicing acceptor creation"


splicing_mutation[
  splicing_mutation$Splicing_Class == "intronic-alternative-5'-splice-site",
  "Splicing_Class"] <- "alternative-5'-splice-site"

splicing_mutation[
  splicing_mutation$Splicing_Class == "intronic-alternative-3'-splice-site",
  "Splicing_Class"] <- "alternative-3'-splice-site"

splicing_mutation[
  splicing_mutation$Splicing_Class == "opposite-side-intron-retention",
  "Splicing_Class"] <- "intron-retention"



splicing_mutation$Mutation_Type <- 
  factor(splicing_mutation$Mutation_Type,
         levels = rev(c("splicing donor disruption",
                    "splicing donor creation",
                    "splicing acceptor disruption",
                    "splicing acceptor creation")))


splicing_mutation_count <- splicing_mutation %>% 
  group_by(Mutation_Type, Splicing_Class) %>% 
  summarize(count = n())

splicing_mutation2 <- splicing_mutation %>% 
  group_by(Mutation_Type, Splicing_Class) %>% 
  summarize(count = n())



ggplot(splicing_mutation2, aes(x = Mutation_Type, y = count, fill = Splicing_Class)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "splicing event count", fill = "splicing class") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))
        ) +
  scale_fill_brewer(palette = "Pastel1")


ggsave("../matome/category_count.png", width = 12, height = 4)


  

library(dplyr)
library(ggplot2)

source("subscript_matome/plot_config.R")

##########
# caterogy summary plot

splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)


splicing_mutation$Mutation_Type <- splicing_mutation$Mutation_Type

splicing_mutation$Is_Indel <- 
  unlist(lapply(strsplit(splicing_mutation$Mutation_Key, ','), 
                function(x) {(nchar(x[3]) > 1 | nchar(x[4]) > 1)}))


splicing_mutation$Mutation_Type2 <- rep("", nrow(splicing_mutation))
splicing_mutation$Mutation_Type3 <- rep("", nrow(splicing_mutation))

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor disruption" &
  splicing_mutation$Is_Canonical == "canonical" &
  splicing_mutation$Is_Indel == FALSE] <- "SNV, canonical donor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, noncanonical donor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor disruption" &
    splicing_mutation$Is_Canonical == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, canonical donor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, noncanonical donor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, donor creation"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, donor creation"


splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, canonical acceptor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, noncanonical acceptor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, canonical acceptor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, noncanonical acceptor disruption"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, acceptor creation"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, acceptor creation"


splicing_mutation$Mutation_Type3[
    splicing_mutation$Mutation_Type2 %in% c("SNV, canonical donor disruption", "SNV, noncanonical donor disruption", 
                                            "Indel, canonical donor disruption", "Indel, noncanonical donor disruption")] <- "Donor disruption"

splicing_mutation$Mutation_Type3[
    splicing_mutation$Mutation_Type2 %in% c("SNV, donor creation", "Indel, donor creation")] <- "Donor creation"

splicing_mutation$Mutation_Type3[
    splicing_mutation$Mutation_Type2 %in% c("SNV, canonical acceptor disruption", "SNV, noncanonical acceptor disruption",
                                            "Indel, canonical acceptor disruption", "Indel, noncanonical acceptor disruption")] <- "Acceptor disruption"

splicing_mutation$Mutation_Type3[
    splicing_mutation$Mutation_Type2 %in% c("SNV, acceptor creation", "Indel, acceptor creation")] <- "Acceptor creation"



splicing_mutation[
  splicing_mutation$Splicing_Class == "intronic-alternative-5'-splice-site",
  "Splicing_Class"] <- "alternative-5'-splice-site"

splicing_mutation[
  splicing_mutation$Splicing_Class == "intronic-alternative-3'-splice-site",
  "Splicing_Class"] <- "alternative-3'-splice-site"

splicing_mutation[
  splicing_mutation$Splicing_Class == "opposite-side-intron-retention",
  "Splicing_Class"] <- "intron-retention"


splicing_mutation$Splicing_Class <-
  factor(splicing_mutation$Splicing_Class,
         levels = c("exon-skip", "alternative-5'-splice-site", 
                    "alternative-3'-splice-site",
                    "intron-retention"),
         labels = c("Exon skip", "Alternative 5' splice site",
                    "Alternative 3' splice site", "Intron retention"))


splicing_mutation$Mutation_Type2 <- 
  factor(splicing_mutation$Mutation_Type2,
         levels = rev(c("SNV, canonical donor disruption",
                        "SNV, noncanonical donor disruption",
                        "SNV, donor creation",
                        "Indel, canonical donor disruption",
                        "Indel, noncanonical donor disruption",
                        "Indel, donor creation",
                        "SNV, canonical acceptor disruption",
                        "SNV, noncanonical acceptor disruption",
                        "SNV, acceptor creation",
                        "Indel, canonical acceptor disruption",
                        "Indel, noncanonical acceptor disruption",
                        "Indel, acceptor creation")))

splicing_mutation$Mutation_Type3 <-
    factor(splicing_mutation$Mutation_Type3,
        levels = rev(c("Donor disruption", "Donor creation", "Acceptor disruption", "Acceptor creation")))


splicing_mutation_count <- splicing_mutation %>% 
  group_by(Mutation_Type2, Splicing_Class) %>% 
  summarize(count = n())

ggplot(splicing_mutation_count, aes(x = Mutation_Type2, y = count, fill = Splicing_Class)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "Splicing event count", fill = "") +
  my_theme() +
  theme(legend.position = "bottom",
  ) +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


ggsave("../matome/category_count.pdf", width = 9, height = 5)



splicing_mutation_count_simple <- splicing_mutation %>%
  group_by(Mutation_Type3, Splicing_Class) %>%
  summarize(count = n())

ggplot(splicing_mutation_count_simple, aes(x = Mutation_Type3, y = count, fill = Splicing_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "Splicing event count", fill = "") +
  my_theme() +
  theme(legend.position = "bottom",
        ) +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave("../matome/category_count_simple.pdf", width = 7, height = 3)


write.table(splicing_mutation %>% 
  select(Cancer_Type, Sample_Name, Mutation_Key, Mutation_Type2) %>% 
  distinct() %>% group_by(Mutation_Type2) %>% summarize(count = n()),
  "../matome/mutation_caterogy_count.txt", quote = FALSE, row.names = FALSE, sep = "\t")



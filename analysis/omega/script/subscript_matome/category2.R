library(dplyr)
library(ggplot2)

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
  splicing_mutation$Is_Indel == FALSE] <- "SNV, donor disruption (C)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, donor disruption (N)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor disruption" &
    splicing_mutation$Is_Canonical == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, donor disruption (C)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, donor disruption (N)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, donor creation"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing donor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, donor creation"


splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, acceptor disruption (C)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, acceptor disruption (N)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, acceptor disruption (C)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, acceptor disruption (N)"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, acceptor creation"

splicing_mutation$Mutation_Type2[
  splicing_mutation$Mutation_Type == "splicing acceptor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, acceptor creation"


splicing_mutation$Mutation_Type3[
    splicing_mutation$Mutation_Type2 %in% c("SNV, donor disruption (C)", "SNV, donor disruption (N)", 
                                            "Indel, donor disruption (C)", "Indel, donor disruption (N)")] <- "Donor disruption"

splicing_mutation$Mutation_Type3[
    splicing_mutation$Mutation_Type2 %in% c("SNV, donor creation", "Indel, donor creation")] <- "Donor creation"

splicing_mutation$Mutation_Type3[
    splicing_mutation$Mutation_Type2 %in% c("SNV, acceptor disruption (C)", "SNV, acceptor disruption (N)",
                                            "Indel, acceptor disruption (C)", "Indel, acceptor disruption (N)")] <- "Acceptor disruption"

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
                    "intron-retention"))


splicing_mutation$Mutation_Type2 <- 
  factor(splicing_mutation$Mutation_Type2,
         levels = rev(c("SNV, donor disruption (C)",
                        "SNV, donor disruption (N)",
                        "SNV, donor creation",
                        "Indel, donor disruption (C)",
                        "Indel, donor disruption (N)",
                        "Indel, donor creation",
                        "SNV, acceptor disruption (C)",
                        "SNV, acceptor disruption (N)",
                        "SNV, acceptor creation",
                        "Indel, acceptor disruption (C)",
                        "Indel, acceptor disruption (N)",
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
  labs(x = "", y = "splicing event count", fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1))
        ) +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

# ggsave("category_count.png", width = 8, height = 4)
ggsave("../matome/category_count.png", width = 12, height = 4)


splicing_mutation_count_simple <- splicing_mutation %>%
  group_by(Mutation_Type3, Splicing_Class) %>%
  summarize(count = n())

ggplot(splicing_mutation_count_simple, aes(x = Mutation_Type3, y = count, fill = Splicing_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "splicing event count", fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.title = element_text(size = rel(1))
        ) +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

# ggsave("category_count.png", width = 8, height = 4)
ggsave("../matome/category_count_simple.png", width = 12, height = 3)


write.table(splicing_mutation %>% 
  select(Cancer_Type, Sample_Name, Mutation_Key, Mutation_Type2) %>% 
  distinct() %>% group_by(Mutation_Type2) %>% summarize(count = n()),
  "../matome/mutation_caterogy_count.txt", quote = FALSE, row.names = FALSE, sep = "\t")


  

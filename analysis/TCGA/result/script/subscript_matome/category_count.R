library(dplyr)
library(ggplot2)

source("../../../conf/plot_config.R")
##########
# caterogy summary plot

splicing_mutation <- read.table("../../output/savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt", header = TRUE, sep = "\t", 
                                quote="", stringsAsFactors = FALSE) %>%
  filter(IR_filtered == "FALSE")

splicing_mutation_proc <- splicing_mutation %>% 
  group_by(Sample_Name, Mutation_Key) %>% 
  summarize(Mutation_Type_Paste = paste(unique(Mutation_Type), collapse = ";"),
            Is_Canonical_Paste = paste(unique(Is_Canonical), collapse = ";")
            )

##########
# processing rule
# canonical > non-canonical
# disruption > creation

select_mut_type <- function(mut_type_str) {
  
  mut_type_vec <- strsplit(mut_type_str, ";")[[1]]
  donor_type_ind <- grep("disruption", mut_type_vec)
  if (length(donor_type_ind) > 0) {
    return(mut_type_vec[min(donor_type_ind)])
  } else {
    return(mut_type_vec[1])
  }
  
}

splicing_mutation_proc$Is_Canonical2 <- splicing_mutation_proc$Is_Canonical_Paste
splicing_mutation_proc$Is_Canonical2[
  splicing_mutation_proc$Is_Canonical2 %in% c("canonical;non-canonical", "non-canonical;canonical")
] <- "canonical"

splicing_mutation_proc$Mutation_Type2 <- unlist(lapply(splicing_mutation_proc$Mutation_Type_Paste, select_mut_type))



##########

splicing_mutation <- splicing_mutation %>% 
  left_join(splicing_mutation_proc, key = c("Sample_Name", "Mutation_Key"))


splicing_mutation$Is_Indel <- 
  unlist(lapply(strsplit(splicing_mutation$Mutation_Key, ','), 
                function(x) {(nchar(x[3]) > 1 | nchar(x[4]) > 1)}))


splicing_mutation$Mutation_Type3 <- rep("", nrow(splicing_mutation))
splicing_mutation$Mutation_Type4 <- rep("", nrow(splicing_mutation))

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing donor disruption" &
  splicing_mutation$Is_Canonical2 == "canonical" &
  splicing_mutation$Is_Indel == FALSE] <- "SNV, canonical donor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing donor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, noncanonical donor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing donor disruption" &
    splicing_mutation$Is_Canonical2 == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, canonical donor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing donor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, noncanonical donor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing donor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, donor creation"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing donor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, donor creation"


splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, canonical acceptor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, noncanonical acceptor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, canonical acceptor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, noncanonical acceptor disruption"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing acceptor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, acceptor creation"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "splicing acceptor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, acceptor creation"


splicing_mutation$Mutation_Type4[
    splicing_mutation$Mutation_Type3 %in% c("SNV, canonical donor disruption", "SNV, noncanonical donor disruption", 
                                            "Indel, canonical donor disruption", "Indel, noncanonical donor disruption")] <- "Donor disruption"

splicing_mutation$Mutation_Type4[
    splicing_mutation$Mutation_Type3 %in% c("SNV, donor creation", "Indel, donor creation")] <- "Donor creation"

splicing_mutation$Mutation_Type4[
    splicing_mutation$Mutation_Type3 %in% c("SNV, canonical acceptor disruption", "SNV, noncanonical acceptor disruption",
                                            "Indel, canonical acceptor disruption", "Indel, noncanonical acceptor disruption")] <- "Acceptor disruption"

splicing_mutation$Mutation_Type4[
    splicing_mutation$Mutation_Type3 %in% c("SNV, acceptor creation", "Indel, acceptor creation")] <- "Acceptor creation"



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
         labels = c("Exon skip", "Alternative 5'-ss",
                    "Alternative 3'-ss", "Intron retention"))


splicing_mutation$Mutation_Type3 <- 
  factor(splicing_mutation$Mutation_Type3,
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

splicing_mutation$Mutation_Type4 <-
    factor(splicing_mutation$Mutation_Type4,
        levels = rev(c("Donor disruption", "Donor creation", "Acceptor disruption", "Acceptor creation")))


splicing_mutation_count <- splicing_mutation %>% 
  group_by(Mutation_Type3, Splicing_Class) %>% 
  summarize(count = n())

ggplot(splicing_mutation_count, aes(x = Mutation_Type3, y = count, fill = Splicing_Class)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "Splicing event count", fill = "") +
  my_theme() +
  theme(legend.position = "bottom",
  ) +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


ggsave("../figure/category_count.pdf", width = 9, height = 5)



splicing_mutation_count_simple <- splicing_mutation %>%
  group_by(Mutation_Type4, Splicing_Class) %>%
  summarize(count = n())

ggplot(splicing_mutation_count_simple, aes(x = Mutation_Type4, y = count, fill = Splicing_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "Splicing event count", fill = "") +
  my_theme() +
  theme(legend.position = "bottom",
        ) +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave("../figure/category_count_simple.pdf", width = 7, height = 3)


write.table(splicing_mutation %>% 
  select(Cancer_Type, Sample_Name, Mutation_Key, Mutation_Type3) %>% 
  distinct() %>% group_by(Mutation_Type3) %>% summarize(count = n()),
  "../temporary/mutation_caterogy_count.txt", quote = FALSE, row.names = FALSE, sep = "\t")



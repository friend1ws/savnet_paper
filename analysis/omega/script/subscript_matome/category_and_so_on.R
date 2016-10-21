library(dplyr)
library(ggplot2)

##########
# caterogy summary plot
splicing_mutation <- read.table("omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)


splicing_mutation[
  intersect(grep("splicing acceptor disruption", splicing_mutation$Motif_Type),
            grep("splicing acceptor creation", splicing_mutation$Motif_Type)), 
  "Motif_Type"] <- "splicing acceptor disruption + creation"

splicing_mutation[
  intersect(grep("splicing donor disruption", splicing_mutation$Motif_Type),
            grep("splicing donor creation", splicing_mutation$Motif_Type)), 
  "Motif_Type"] <- "splicing donor disruption + creation"


tind <- which(splicing_mutation$Motif_Type == "splicing acceptor disruption;splicing acceptor disruption")
is_canonical_ind <- unlist(lapply(strsplit(splicing_mutation[tind, "Is_Canonical"], ";"), 
                                  function(x) {sum(x == "canonical") > 0}))

splicing_mutation[tind[is_canonical_ind], "Is_Canonical"] <- "canonical"
splicing_mutation[tind[!is_canonical_ind], "Is_Canonical"] <- "non-canonical"
splicing_mutation[tind, "Motif_Type"] <- "splicing acceptor disruption"


tind <- which(splicing_mutation$Motif_Type == "splicing acceptor creation;splicing acceptor creation")
is_canonical_ind <- unlist(lapply(strsplit(splicing_mutation[tind, "Is_Canonical"], ";"), 
                                  function(x) {sum(x == "canonical") > 0}))

splicing_mutation[tind[is_canonical_ind], "Is_Canonical"] <- "canonical"
splicing_mutation[tind[!is_canonical_ind], "Is_Canonical"] <- "non-canonical"
splicing_mutation[tind, "Motif_Type"] <- "splicing acceptor creation"



splicing_mutation$Mutation_Type <- splicing_mutation$Motif_Type


splicing_mutation[
  splicing_mutation$Splicing_Type == "intronic-alternative-5'-splice-site",
  "Splicing_Type"] <- "alternative-5'-splice-site"

splicing_mutation[
  splicing_mutation$Splicing_Type == "intronic-alternative-3'-splice-site",
  "Splicing_Type"] <- "alternative-3'-splice-site"


cn_ind <- splicing_mutation$Motif_Type %in% c("splicing donor disruption", "splicing donor creation",
                                              "splicing acceptor disruption", "splicing acceptor creation")

splicing_mutation[cn_ind & splicing_mutation$Is_Canonical == "canonical", "Mutation_Type"] <- 
  paste(splicing_mutation[cn_ind & splicing_mutation$Is_Canonical == "canonical", "Mutation_Type"] , "(C)")

splicing_mutation[cn_ind & splicing_mutation$Is_Canonical == "non-canonical", "Mutation_Type"] <- 
  paste(splicing_mutation[cn_ind & splicing_mutation$Is_Canonical == "non-canonical", "Mutation_Type"] , "(N)")


splicing_mutation$Mutation_Type <- 
  factor(splicing_mutation$Mutation_Type,
         levels = c("splicing donor disruption (C)", "splicing donor disruption (N)",
                    "splicing donor creation (C)", "splicing donor creation (N)",
                    "splicing acceptor disruption (C)", "splicing acceptor disruption (N)",
                    "splicing acceptor creation (C)", "splicing acceptor creation (N)",
                    "splicing donor disruption + creation", "splicing acceptor disruption + creation"))


splicing_mutation2 <- splicing_mutation %>% 
  group_by(Mutation_Type, Splicing_Type) %>% 
  summarize(count = n())



ggplot(splicing_mutation2, aes(x = Mutation_Type, y = count, fill = Splicing_Type)) + 
  geom_bar(stat = "identity") +
  coord_flip()


splicing_mutation3 <- splicing_mutation %>% 
  group_by(Cancer_Type, Mutation_Type) %>% 
  summarize(count = n())


ggplot(splicing_mutation3, aes(x = Cancer_Type, y = count, fill = Mutation_Type)) + 
  geom_bar(stat = "identity") +
  coord_flip()

  
##########
splicing_mutation <- read.table("omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

multiple_effect <- splicing_mutation %>% 
  group_by(Cancer_Type, Sample_Name, Mutation_Key) %>% 
  summarize(splice_count = n()) %>%
  group_by(splice_count) %>% summarize(count = n())


ggplot(multiple_effect, aes(x = splice_count, y = count)) + 
  geom_bar(stat = "identity") + 
  scale_y_log10() +
  theme_minimal()






##########
# cancer gene summary
splicing_mutation <- read.table("omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

cancer_splicing_mutation <- splicing_mutation %>% 
  filter(Is_Cancer_Gene == "TRUE") %>% 
  select(Cancer_Type, Sample_Name, Gene_Symbol) %>%
  distinct() %>%
  group_by(Gene_Symbol) %>%
  summarize(count = n()) %>% arrange(desc(count))


##########
# hot spot mutation
splicing_mutation <- read.table("omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

hotspot_splicing_mutation <- splicing_mutation %>% 
  select(Cancer_Type, Sample_Name, Gene_Symbol, Mutation_Key, Is_Canonical) %>%
  distinct() 

hs_mut_count <- hotspot_splicing_mutation %>% 
  group_by(Mutation_Key) %>%
  summarize(mutation_count = n()) 

hotspot_splicing_mutation <- left_join(hotspot_splicing_mutation, hs_mut_count, by = "Mutation_Key") %>%
  filter(mutation_count >= 3)

hotsplot_splicing_mutation_count <- hotspot_splicing_mutation %>% 
  select(Gene_Symbol, Mutation_Key, Is_Canonical, mutation_count) %>% 
  distinct() %>% arrange(desc(mutation_count))

##########
# hot spot splicing
splicing_mutation <- read.table("omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

hotspot_splicing_key <- splicing_mutation %>% 
  select(Cancer_Type, Sample_Name, Gene_Symbol, Splicing_Key, Splicing_Type, Is_Inframe) %>%
  distinct() 

hs_splicing_count <- hotspot_splicing_key %>% 
  group_by(Splicing_Key) %>%
  summarize(splicing_count = n()) 

hotspot_splicing_key <- left_join(hotspot_splicing_key, hs_splicing_count, by = "Splicing_Key") %>%
  filter(splicing_count >= 3)

hotspot_splicing_key_count <- hotspot_splicing_key %>% 
  select(Gene_Symbol, Splicing_Key, Splicing_Type, Is_Inframe, splicing_count) %>% 
  distinct() %>% arrange(desc(splicing_count))


##########
# inframe info

inframe_info <- read.table("inframe_summary.txt", sep = "\t", header = TRUE)






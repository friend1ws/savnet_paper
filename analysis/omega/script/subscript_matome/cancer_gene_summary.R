library(dplyr)
library(ggplot2)

##########
# gene summary
splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

# cancer_gene
cancer_gene <- splicing_mutation %>% 
  filter(Is_Cancer_Gene == "TRUE") %>% 
  select(Gene_Symbol) %>% 
  distinct()

splicing_mutation_count <- splicing_mutation %>% 
  # filter(Is_Cancer_Gene == "TRUE") %>% 
  select(Cancer_Type, Sample_Name, Gene_Symbol) %>%
  distinct() %>%
  group_by(Cancer_Type, Gene_Symbol) %>%
  summarize(count = n()) %>% arrange(desc(count))

splicing_mutation_count_total <- splicing_mutation_count %>%
  group_by(Gene_Symbol) %>%
  summarize(total_count = sum(count)) %>%
  filter(total_count >= 10) %>% 
  arrange(total_count)


splicing_mutation_count_proc <- splicing_mutation_count %>% 
  filter(Gene_Symbol %in% splicing_mutation_count_total$Gene_Symbol)

splicing_mutation_count_proc$Gene_Symbol2 <-
  factor(splicing_mutation_count_proc$Gene_Symbol,
         levels = splicing_mutation_count_total$Gene_Symbol)

gene_colour <- rep("grey30", nrow(splicing_mutation_count_proc))
gene_colour[levels(splicing_mutation_count_proc$Gene_Symbol2) %in% cancer_gene$Gene_Symbol] <- "red"

ggplot(splicing_mutation_count_proc,
       aes(x = Cancer_Type, y = Gene_Symbol2, size = count)) + 
  geom_point() +
  theme_minimal() +
  labs(x = "Cancer Type", y = "Gene") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2), colour =  gene_colour),
        axis.title = element_text(size = rel(1.2)))


ggsave("../matome/gene_cancertype_summary.png", width = 12, height = 8, units = "in")


##########
# ratio

splicing_mutation <- read.table("../matome/omega.genomon_splicing_mutation.result.txt", header = TRUE, sep = "\t", as.is=TRUE, quote="", stringsAsFactors = FALSE)

splicing_mutation <- splicing_mutation %>% 
  filter(Mutation_Type %in% c("splicing acceptor disruption",
                             "splicing donor disruption",
                             "splicing donor creation",
                             "splicing acceptor creation"))

# splicing_mutation$Mutation_Type <- rep("", nrow(splicing_mutation))

splicing_mutation[
  grep("creation", splicing_mutation$Mutation_Type), 
  "Mutation_Type"] <- "motif creation"

splicing_mutation[
  intersect(grep("disruption", splicing_mutation$Mutation_Type),
            which(splicing_mutation$Is_Canonical == "canonical")),
  "Mutation_Type"] <- "motif disruption (C)"

splicing_mutation[
  intersect(grep("disruption", splicing_mutation$Mutation_Type),
            which(splicing_mutation$Is_Canonical == "non-canonical")),
  "Mutation_Type"] <- "motif disruption (N)"



cancer_splicing_mutation_count <- splicing_mutation %>% 
  filter(Is_Cancer_Gene == "TRUE") %>% 
  select(Cancer_Type, Sample_Name, Gene_Symbol, Mutation_Type) %>%
  distinct() %>%
  group_by(Gene_Symbol, Mutation_Type) %>%
  summarize(count = n()) 

cancer_splicing_mutation_count_total <- cancer_splicing_mutation_count %>%
  group_by(Gene_Symbol) %>%
  summarize(total_count = sum(count)) %>%
  filter(total_count >= 8) %>% 
  arrange(total_count)


cancer_splicing_mutation_count_proc <- cancer_splicing_mutation_count %>% 
  filter(Gene_Symbol %in% cancer_splicing_mutation_count_total$Gene_Symbol)

cancer_splicing_mutation_count_proc$Gene_Symbol2 <-
  factor(cancer_splicing_mutation_count_proc$Gene_Symbol,
         levels = cancer_splicing_mutation_count_total$Gene_Symbol)

cancer_splicing_mutation_count_proc$Mutation_Type2 <-
  factor(cancer_splicing_mutation_count_proc$Mutation_Type,
         levels = c("motif disruption (C)",
                    "motif disruption (N)",
                    "motif creation"))


ggplot(cancer_splicing_mutation_count_proc %>% arrange(Mutation_Type2), 
       aes(x = Gene_Symbol2, y = count, fill = Mutation_Type2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "gene", y = "#mutation", fill = "") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = rel(1.0)),
        axis.title = element_text(size = rel(1.0))) +
  scale_fill_brewer(palette = "Pastel2")

ggsave("../matome/gene_mutationtype_summary.png", width = 10, height = 6, units = "in")

##########
##########
# inframe info

inframe_info <- read.table("../matome/gene.inframe_summary.txt", sep = "\t", header = TRUE)



inframe_count <- inframe_info %>% 
  filter(Is_Cancer_Gene == "TRUE") %>% 
  select(Cancer_Type, Sample_Name, Gene_Symbol, Mutation_Key, Is_Inframe) %>%
  distinct() %>%
  group_by(Gene_Symbol, Is_Inframe) %>%
  summarize(count = n())

inframe_count_total <- inframe_count %>%
  group_by(Gene_Symbol) %>% 
  summarize(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>% 
  filter(total_count >= 8) %>%
  arrange(total_count)


inframe_count_proc <- inframe_count %>% 
  filter(Gene_Symbol %in% inframe_count_total$Gene_Symbol)

inframe_count_proc$Gene_Symbol2 <-
  factor(inframe_count_proc$Gene_Symbol,
         levels = inframe_count_total$Gene_Symbol)


ggplot(inframe_count_proc, 
       aes(x = Gene_Symbol2, y = count, fill = Is_Inframe)) + 
  geom_bar(stat = "identity") +
  labs(x = "gene", y = "#mutation", fill = "inframe") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2")

ggsave("../matome/inframe_summary.png", width = 7, height = 4)




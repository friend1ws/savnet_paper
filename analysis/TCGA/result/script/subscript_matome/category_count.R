library(dplyr)
library(ggplot2)
library(cowplot)

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


splicing_mutation_proc$Mutation_Type2 <-
  factor(splicing_mutation_proc$Mutation_Type2,
         levels = rev(c("splicing donor disruption", "splicing donor creation",
                        "splicing acceptor disruption", "splicing acceptor creation")),
         labels = rev(c("Donor disruption", "Donor creation", 
                        "Acceptor disruption", "Acceptor creation")))


##########


splicing_mutation <- splicing_mutation %>% 
  left_join(splicing_mutation_proc, key = c("Sample_Name", "Mutation_Key"))


splicing_mutation$Is_Indel <- 
  unlist(lapply(strsplit(splicing_mutation$Mutation_Key, ','), 
                function(x) {(nchar(x[3]) > 1 | nchar(x[4]) > 1)}))


splicing_mutation$Mutation_Type3 <- rep("", nrow(splicing_mutation))
# splicing_mutation$Mutation_Type4 <- rep("", nrow(splicing_mutation))

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Donor disruption" &
  splicing_mutation$Is_Canonical2 == "canonical" &
  splicing_mutation$Is_Indel == FALSE] <- "SNV, canonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Donor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, noncanonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Donor disruption" &
    splicing_mutation$Is_Canonical2 == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, canonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Donor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, noncanonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Donor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Donor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel"


splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, canonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV, noncanonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, canonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Acceptor disruption" &
    splicing_mutation$Is_Canonical2 == "non-canonical" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel, noncanonical"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Acceptor creation" &
    splicing_mutation$Is_Indel == FALSE] <- "SNV"

splicing_mutation$Mutation_Type3[
  splicing_mutation$Mutation_Type2 == "Acceptor creation" &
    splicing_mutation$Is_Indel == TRUE] <- "Indel"


# splicing_mutation$Mutation_Type4[
#     splicing_mutation$Mutation_Type3 %in% c("SNV, canonical donor disruption", "SNV, noncanonical donor disruption", 
#                                             "Indel, canonical donor disruption", "Indel, noncanonical donor disruption")] <- "Donor disruption"

# splicing_mutation$Mutation_Type4[
#     splicing_mutation$Mutation_Type3 %in% c("SNV, donor creation", "Indel, donor creation")] <- "Donor creation"

# splicing_mutation$Mutation_Type4[
#     splicing_mutation$Mutation_Type3 %in% c("SNV, canonical acceptor disruption", "SNV, noncanonical acceptor disruption",
#                                             "Indel, canonical acceptor disruption", "Indel, noncanonical acceptor disruption")] <- "Acceptor disruption"

# splicing_mutation$Mutation_Type4[
#     splicing_mutation$Mutation_Type3 %in% c("SNV, acceptor creation", "Indel, acceptor creation")] <- "Acceptor creation"




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
         labels = c("Exon skipping", "Alternative 5'SS",
                    "Alternative 3'SS", "Intron retention"))


mut_type3_order <- rev(c("SNV, canonical", "SNV, noncanonical", "SNV", "Indel, canonical", "Indel, noncanonical", "Indel"))

splicing_mutation$Mutation_Type3 <- 
  factor(splicing_mutation$Mutation_Type3,
         levels = mut_type3_order)

         

# splicing_mutation$Mutation_Type4 <-
#     factor(splicing_mutation$Mutation_Type4,
#         levels = rev(c("Donor disruption", "Donor creation", "Acceptor disruption", "Acceptor creation")))


# splicing_mutation_count <- splicing_mutation %>% 
#   group_by(Mutation_Type3, Splicing_Class) %>% 
#  summarize(count = n())

splicing_mutation_type_count2 <- splicing_mutation %>%
    filter(Mutation_Type2 == "Donor disruption") %>%
    group_by(Mutation_Type3, Splicing_Class) %>%
    summarize(count = n())

splicing_mutation_type_count3 <- splicing_mutation %>%
    filter(Mutation_Type2 == "Donor disruption") %>%
    select(Sample_Name, Mutation_Key, Mutation_Type3) %>%
    distinct() %>%
    group_by(Mutation_Type3) %>%
    summarize(count = n())

splicing_mutation_type_count2$Mutation_Type32 <- factor(splicing_mutation_type_count2$Mutation_Type3,
    levels = mut_type3_order,
    labels = unlist(lapply(mut_type3_order, 
                           function(x) {paste(x, " (", splicing_mutation_type_count3$count[splicing_mutation_type_count3$Mutation_Type3 == x], ")", sep = "")})))


print(splicing_mutation_type_count3)


# p_dd <- ggplot(splicing_mutation %>% 
#          filter(Mutation_Type2 == "Donor disruption") %>% 
#          group_by(Mutation_Type3, Splicing_Class) %>% 
#          summarize(count = n()), 
#        aes(x = Mutation_Type3, y = count, fill = Splicing_Class)) + 
p_dd <-ggplot(splicing_mutation_type_count2, aes(Mutation_Type32, y = count, fill = Splicing_Class)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle(paste("Donor disruption (", sum(splicing_mutation_type_count3$count), ")", sep ="")) +
  labs(x = "", y = "", fill = "") +
  my_theme() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5300)) +
  guides(fill = FALSE)
  # guides(fill=guide_legend(nrow=1,byrow=TRUE))


splicing_mutation_type_count2 <- splicing_mutation %>%
    filter(Mutation_Type2 == "Donor creation") %>%
    group_by(Mutation_Type3, Splicing_Class) %>%
    summarize(count = n())

splicing_mutation_type_count3 <- splicing_mutation %>%
    filter(Mutation_Type2 == "Donor creation") %>%
    select(Sample_Name, Mutation_Key, Mutation_Type3) %>%
    distinct() %>%
    group_by(Mutation_Type3) %>%
    summarize(count = n())

splicing_mutation_type_count2$Mutation_Type32 <- factor(splicing_mutation_type_count2$Mutation_Type3,
    levels = mut_type3_order,
    labels = unlist(lapply(mut_type3_order,
                           function(x) {paste(x, " (", splicing_mutation_type_count3$count[splicing_mutation_type_count3$Mutation_Type3 == x], ")", sep = "")})))


p_dc <- ggplot(splicing_mutation_type_count2, aes(Mutation_Type32, y = count, fill = Splicing_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle(paste("Donor creation (", sum(splicing_mutation_type_count3$count), ")", sep ="")) +
  labs(x = "", y = "Splicing event count", fill = "") +
  my_theme() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 5300)) +
  guides(fill = FALSE)
# guides(fill=guide_legend(nrow=1,byrow=TRUE))


splicing_mutation_type_count2 <- splicing_mutation %>% 
    filter(Mutation_Type2 == "Acceptor disruption") %>%
    group_by(Mutation_Type3, Splicing_Class) %>% 
    summarize(count = n())

splicing_mutation_type_count3 <- splicing_mutation %>% 
    filter(Mutation_Type2 == "Acceptor disruption") %>%
    select(Sample_Name, Mutation_Key, Mutation_Type3) %>% 
    distinct() %>% 
    group_by(Mutation_Type3) %>% 
    summarize(count = n())

splicing_mutation_type_count2$Mutation_Type32 <- factor(splicing_mutation_type_count2$Mutation_Type3,
    levels = mut_type3_order,
    labels = unlist(lapply(mut_type3_order,
                           function(x) {paste(x, " (", splicing_mutation_type_count3$count[splicing_mutation_type_count3$Mutation_Type3 == x], ")", sep = "")})))


p_ad <- ggplot(splicing_mutation_type_count2, aes(Mutation_Type32, y = count, fill = Splicing_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle(paste("Acceptor disruption (", sum(splicing_mutation_type_count3$count), ")", sep ="")) +
  labs(x = "", y = "", fill = "") +
  my_theme() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5300)) +
  guides(fill = FALSE)
# guides(fill=guide_legend(nrow=1,byrow=TRUE))


splicing_mutation_type_count2 <- splicing_mutation %>%
    filter(Mutation_Type2 == "Acceptor creation") %>%
    group_by(Mutation_Type3, Splicing_Class) %>%
    summarize(count = n())

splicing_mutation_type_count3 <- splicing_mutation %>%
    filter(Mutation_Type2 == "Acceptor creation") %>%
    select(Sample_Name, Mutation_Key, Mutation_Type3) %>%
    distinct() %>%
    group_by(Mutation_Type3) %>%
    summarize(count = n())

splicing_mutation_type_count2$Mutation_Type32 <- factor(splicing_mutation_type_count2$Mutation_Type3,
    levels = mut_type3_order,
    labels = unlist(lapply(mut_type3_order,
                           function(x) {paste(x, " (", splicing_mutation_type_count3$count[splicing_mutation_type_count3$Mutation_Type3 == x], ")", sep = "")})))


p_ac <- ggplot(splicing_mutation_type_count2, aes(Mutation_Type32, y = count, fill = Splicing_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle(paste("Acceptor creation (", sum(splicing_mutation_type_count3$count), ")", sep ="")) +
  labs(x = "", y = "Splicing event count", fill = "") +
  my_theme() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 5300)) +
  guides(fill = FALSE)
# guides(fill=guide_legend(nrow=1,byrow=TRUE))


p_dummy_for_legend <- 
  ggplot(splicing_mutation %>% 
           filter(Mutation_Type2 == "Donor disruption") %>% 
           group_by(Mutation_Type3, Splicing_Class) %>% 
           summarize(count = n()), 
         aes(x = Mutation_Type3, y = count, fill = Splicing_Class)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("Donor disruption") +
  labs(x = "", y = "", fill = "") +
  my_theme() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5300))


p_dd_ad_dc_ac <- plot_grid(p_dd, p_ad, p_dc, p_ac, ncol = 2, align = "hv", rel_heights = c(1, 0.65))

plot_grid(p_dd_ad_dc_ac, g_legend(p_dummy_for_legend), ncol = 1, rel_heights = c(1, 0.1))



ggsave("../figure/category_count.tiff", width = 18, height = 8, dpi = 600, units = "cm")



splicing_mutation_count_simple <- splicing_mutation %>%
  group_by(Mutation_Type2, Splicing_Class) %>%
  summarize(count = n())

ggplot(splicing_mutation_count_simple, aes(x = Mutation_Type2, y = count, fill = Splicing_Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "", y = "Splicing event count", fill = "") +
  my_theme() +
  theme(legend.position = "bottom",
        plot.margin = unit(c(5.5, 7.5, 5.5, 5.5), "points")) +
  scale_fill_manual(values = splicing_class_colour) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

ggsave("../figure/category_count_simple.tiff", width = 10, height = 4, dpi = 600, units = "cm")


write.table(splicing_mutation %>% 
              mutate(Mutation_Type4 = paste(Mutation_Type2, Mutation_Type3, sep =", ")) %>% 
  select(Cancer_Type, Sample_Name, Mutation_Key, Mutation_Type4) %>% 
  distinct() %>% group_by(Mutation_Type4) %>% summarize(count = n()),
  "../temporary/mutation_caterogy_count.txt", quote = FALSE, row.names = FALSE, sep = "\t")



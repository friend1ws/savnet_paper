library(dplyr)
library(spliceSites)

mes<-load.maxEnt()
hb<-load.hbond()

get_pos_df <- function(motif_pos) {
  
  motif_pos_sp1 <- strsplit(motif_pos, ":")
  tchr <- paste("chr", unlist(lapply(motif_pos_sp1, '[', 1)), sep = "")

  motif_pos_sp2 <- strsplit(unlist(lapply(motif_pos_sp1, '[', 2)), ",")
  tstrand <- unlist(lapply(motif_pos_sp2, '[', 2))

  motif_pos_sp3 <- strsplit(unlist(lapply(motif_pos_sp2, '[', 1)), "-")
  tstart <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 1)))
  tend <- as.numeric(unlist(lapply(motif_pos_sp3, '[', 2)))

  return(data.frame(chr = tchr, start = tstart, end = tend, strand = tstrand))
  
}


args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]


snv_info <- read.table(input_file, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE)

canonical_count_dd <- snv_info %>% 
  filter(Mutation_Type == "splicing donor disruption", Is_Canonical == "canonical") %>% 
  group_by(Rel_Pos) %>% summarize(count = n())

if (length(canonical_count_dd$Rel_Pos) != 2) 
  stop("canonical info of donor splicing motif is inconsistent!")


canonical_count_ad <- snv_info %>% 
  filter(Mutation_Type == "splicing acceptor disruption", Is_Canonical == "canonical") %>% 
  group_by(Rel_Pos) %>% summarize(count = n())

if (length(canonical_count_ad$Rel_Pos) != 2) 
  stop("canonical info of acceptor splicing motif is inconsistent!")



##########
# donor maxEnt score
snv_info_d <- snv_info %>% 
  filter(Mutation_Type %in% c("splicing donor disruption", "splicing donor creation"))

pos_df_d <- get_pos_df(snv_info_d$Motif_Pos)

motif_len_d <- unique(pos_df_d$end - pos_df_d$start + 1)
if (length(motif_len_d) != 1) stop("donor motif size is inconsistent!")


exon_size_d <- min(canonical_count_dd$Rel_Pos) - 1
intron_size_d <- motif_len_d - exon_size_d

ent_start_d <- rep(0, length(pos_df_d$start))
ent_start_d[pos_df_d$strand == "+"] <- pos_df_d$start[pos_df_d$strand == "+"] + exon_size_d - 3
ent_start_d[pos_df_d$strand == "-"] <- pos_df_d$start[pos_df_d$strand == "-"] + intron_size_d - 8

ent_end_d <- rep(0, length(pos_df_d$start))
ent_end_d[pos_df_d$strand == "+"] <- pos_df_d$end[pos_df_d$strand == "+"] - intron_size_d + 8
ent_end_d[pos_df_d$strand == "-"] <- pos_df_d$end[pos_df_d$strand == "-"]- exon_size_d + 3

gr <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = pos_df_d$chr, start = ent_start_d, end = ent_end_d, 
             strand = pos_df_d$strand), ignore.strand = FALSE)

donor_seq_wt <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr)

donor_mes_wt <- unlist(lapply(as.character(donor_seq_wt), function(x) {score5(mes, x, pos = 3)}))

hb_score_wt <- hbond(hb, as.character(donor_seq_wt), 3)


mut_pos_d <- snv_info_d$Rel_Pos + 3 - exon_size_d
mut_alt_d <- snv_info_d$Alt_Base

donor_seq_mt <- unlist(lapply(1:length(donor_seq_wt), 
                              function(x) {
                                paste(substring(donor_seq_wt[x], 1, mut_pos_d[x] - 1), 
                                      mut_alt_d[x], 
                                      substring(donor_seq_wt[x], mut_pos_d[x] + 1, 11), 
                                      sep ="")
                              }))

donor_mes_mt <- unlist(lapply(as.character(donor_seq_mt), function(x) {score5(mes, x, pos = 3)}))
hb_score_mt <- hbond(hb, as.character(donor_seq_mt), 3)
##########

##########
# acceptor masEntScore

snv_info_a <- snv_info %>% 
  filter(Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation"))

pos_df_a <- get_pos_df(snv_info_a$Motif_Pos)

motif_len_a <- unique(pos_df_a$end - pos_df_a$start + 1)
if (length(motif_len_a) != 1) stop("donor motif size is inconsistent!")


intron_size_a <- max(canonical_count_ad$Rel_Pos) 
exon_size_a <- motif_len_a - intron_size_a

ent_start_a <- rep(0, length(pos_df_a$start))
ent_start_a[pos_df_a$strand == "+"] <- pos_df_a$start[pos_df_a$strand == "+"] + intron_size_a - 20 
ent_start_a[pos_df_a$strand == "-"] <- pos_df_a$start[pos_df_a$strand == "-"] + exon_size_a - 3

ent_end_a <- rep(0, length(pos_df_a$start))
ent_end_a[pos_df_a$strand == "+"] <- pos_df_a$end[pos_df_a$strand == "+"] - exon_size_a + 3
ent_end_a[pos_df_a$strand == "-"] <- pos_df_a$end[pos_df_a$strand == "-"]- intron_size_a + 20 

gr <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(chr = pos_df_a$chr, start = ent_start_a, end = ent_end_a, 
             strand = pos_df_a$strand), ignore.strand = FALSE)

acceptor_seq_wt <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr)

acceptor_mes_wt <- unlist(lapply(as.character(acceptor_seq_wt), function(x) {score3(mes, x, pos = 20)}))


mut_pos_a <- snv_info_a$Rel_Pos + 20 - intron_size_a
mut_alt_a <- snv_info_a$Alt_Base

acceptor_seq_mt <- unlist(lapply(1:length(acceptor_seq_wt), 
                              function(x) {
                                paste(substring(acceptor_seq_wt[x], 1, mut_pos_a[x] - 1), 
                                      mut_alt_a[x], 
                                      substring(acceptor_seq_wt[x], mut_pos_a[x] + 1, 23), 
                                      sep ="")
                              }))

acceptor_mes_mt <- unlist(lapply(as.character(acceptor_seq_mt), function(x) {score3(mes, x, pos = 20)}))
##########


mes_wt <- rep(0, nrow(snv_info))
mes_wt[snv_info$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation")] <- donor_mes_wt
mes_wt[snv_info$Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation")] <- acceptor_mes_wt

mes_mt <- rep(0, nrow(snv_info))
mes_mt[snv_info$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation")] <- donor_mes_mt
mes_mt[snv_info$Mutation_Type %in% c("splicing acceptor disruption", "splicing acceptor creation")] <- acceptor_mes_mt

hb_wt <- rep(NA, nrow(snv_info))
hb_wt[snv_info$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation")] <- hb_score_wt

hb_mt <- rep(NA, nrow(snv_info))
hb_mt[snv_info$Mutation_Type %in% c("splicing donor disruption", "splicing donor creation")] <- hb_score_mt

write.table(cbind(snv_info, data.frame(mes_wt = mes_wt, mes_mt = mes_mt, hbond_wt = hb_wt, hbond_mt = hb_mt)), 
            file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")




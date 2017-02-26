#! /usr/local/package/r/3.2.5/bin/Rscript

library(dplyr)

gsm_input <- commandArgs(trailingOnly=TRUE)[1]
gsm_prefix <- commandArgs(trailingOnly=TRUE)[2]
output_prefix <- commandArgs(trailingOnly=TRUE)[3]
ctype <- commandArgs(trailingOnly=TRUE)[4]
sf_gene <- commandArgs(trailingOnly=TRUE)[5]
sf_sample_list <- commandArgs(trailingOnly=TRUE)[6]


mut_list <- read.table(gsm_input, sep ="\t", header = TRUE, stringsAsFactors = FALSE) 
mut_num <- nrow(mut_list)

weight <- mut_list[,2] / mean(mut_list[,2])

SF_mut <- read.table(sf_sample_list, sep = "\t", stringsAsFactors = FALSE)

sf_mut_sample <- SF_mut %>% filter(V2 == ctype) %>% filter(V9 == sf_gene)
sf_mut_id <- which(mut_list[,1] %in% sf_mut_sample$V1)

print(sf_mut_id)
print(sf_mut_sample)
print(setdiff(1:mut_num, sf_mut_id))

asp_count <- read.table(paste(gsm_prefix, ".splicing.associate.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
 

sp_vec_neg <- rep("", nrow(asp_count))
sp_vec_pos <- rep("", nrow(asp_count))
Minus_Log_PV <- rep(0, nrow(asp_count))
Effect_Size <- rep(0, nrow(asp_count))
for (i in 1:nrow(asp_count)) {

  sp_vec <- as.numeric(unlist(strsplit(asp_count$Read_Counts[i], ',')))
  sp_vec <- sp_vec / weight

  T <- t.test(sp_vec[setdiff(1:mut_num, sf_mut_id)], sp_vec[sf_mut_id], alternative = "less")
  Minus_Log_PV[i] <- - log10(T$p.value)
  Effect_Size[i] <- mean(sp_vec[sf_mut_id]) / (mean(sp_vec[setdiff(1:mut_num, sf_mut_id)]) + 0)
  sp_vec_neg[i] <- paste(round(sp_vec[setdiff(1:mut_num, sf_mut_id)], 2), collapse = ",")
  sp_vec_pos[i] <- paste(round(sp_vec[sf_mut_id], 2), collapse = ",")

}

write.table(cbind(asp_count[,c(1:4,6:9)], sp_vec_neg, sp_vec_pos, Minus_Log_PV, Effect_Size), paste(output_prefix, "_", sf_gene, "_SJ_IR_filt.txt", sep = ""), sep ="\t", row.names = FALSE, quote = FALSE)

write.table(cbind(asp_count[,c(1:4,6:9)], sp_vec_neg, sp_vec_pos, Minus_Log_PV, Effect_Size)[Minus_Log_PV >= 2,], paste(output_prefix, "_", sf_gene, "_SJ_IR_filt.sig.txt", sep = ""), sep ="\t", row.names = FALSE, quote = FALSE)


sj_count <- read.table(paste(gsm_prefix, ".SJ_merged.annot.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "") %>% 
    filter(Splicing_Class %in% c("alternative-3'-splice-site", "alternative-5'-splice-site", "exon-skip", 
                                 "intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site"))


sp_vec_neg <- rep("", nrow(sj_count))
sp_vec_pos <- rep("", nrow(sj_count))
Minus_Log_PV <- rep(0, nrow(sj_count))
Effect_Size <- rep(0, nrow(sj_count))
for (i in 1:nrow(sj_count)) {
  
  sp_vec <- as.numeric(unlist(strsplit(sj_count$SJ_4[i], ',')))
  sp_vec <- sp_vec / weight 
  
  T <- t.test(sp_vec[setdiff(1:mut_num, sf_mut_id)], sp_vec[sf_mut_id], alternative = "less")
  Minus_Log_PV[i] <- - log10(T$p.value)
  Effect_Size[i] <- mean(sp_vec[sf_mut_id]) / (mean(sp_vec[setdiff(1:mut_num, sf_mut_id)]) + 0)
  sp_vec_neg[i] <- paste(round(sp_vec[setdiff(1:mut_num, sf_mut_id)], 2), collapse = ",")
  sp_vec_pos[i] <- paste(round(sp_vec[sf_mut_id], 2), collapse = ",")

}

write.table(cbind(sj_count[,c(1:3,5:14)], sp_vec_neg, sp_vec_pos, Minus_Log_PV, Effect_Size), paste(output_prefix, "_", sf_gene, "_SJ.txt", sep = ""), sep ="\t", row.names = FALSE, quote = FALSE)

write.table(cbind(sj_count[,c(1:3,5:14)], sp_vec_neg, sp_vec_pos, Minus_Log_PV, Effect_Size)[Minus_Log_PV >= 2,], paste(output_prefix, "_", sf_gene, "_SJ.sig.txt", sep = ""), sep ="\t", row.names = FALSE, quote = FALSE)


ir_count <- read.table(paste(gsm_prefix, ".IR_merged.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

sp_vec_neg <- rep("", nrow(ir_count))
sp_vec_pos <- rep("", nrow(ir_count))
Minus_Log_PV <- rep(0, nrow(ir_count))
Effect_Size <- rep(0, nrow(ir_count))
for (i in 1:nrow(ir_count)) {

  sp_vec <- as.numeric(unlist(strsplit(ir_count$Read_Count_Vector[i], ',')))
  sp_vec <- sp_vec / weight

  T <- t.test(sp_vec[setdiff(1:mut_num, sf_mut_id)], sp_vec[sf_mut_id], alternative = "less")
  Minus_Log_PV[i] <- - log10(T$p.value)
  Effect_Size[i] <- mean(sp_vec[sf_mut_id]) / (mean(sp_vec[setdiff(1:mut_num, sf_mut_id)]) + 0)
  sp_vec_neg[i] <- paste(round(sp_vec[setdiff(1:mut_num, sf_mut_id)], 2), collapse = ",")
  sp_vec_pos[i] <- paste(round(sp_vec[sf_mut_id], 2), collapse = ",")

}

write.table(cbind(ir_count[,c(1:8)], sp_vec_neg, sp_vec_pos, Minus_Log_PV, Effect_Size), paste(output_prefix, "_", sf_gene, "_IR.txt", sep = ""), sep ="\t", row.names = FALSE, quote = FALSE)

write.table(cbind(ir_count[,c(1:8)], sp_vec_neg, sp_vec_pos, Minus_Log_PV, Effect_Size)[Minus_Log_PV >= 2.0,], paste(output_prefix, "_", sf_gene, "_IR.sig.txt", sep = ""), sep ="\t", row.names = FALSE, quote = FALSE)


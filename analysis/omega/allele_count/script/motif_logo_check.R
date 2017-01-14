library(dplyr)
library(seqLogo)


splicing_mutation <- read.table("../output/omega.splicing_mutation.info.txt", sep = "\t", header = TRUE, quote = "",
                                stringsAsFactors = FALSE) %>%
  filter(Ref_Mut != "-" & Alt_Mut != "-") %>% filter(FKPM >= 10.0)


if (!file.exists("../output/motif_diff")) {
  dir.create("../output/motif_diff")
}
if (!file.exists("../output/motif_diff_ns")) {
  dir.create("../output/motif_diff_ns")
}


for(pos in 1:9) {
  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GenomonSplicingMutation == "---") 
  
  base_count_mat <- matrix(0, 9, 4)
  for(n in 1:9) {
    base_count <- table(substring(motif_sub$Motif_Seq, n + 0, n + 0))
    base_count_mat[n, 1] <- base_count["A"]
    base_count_mat[n, 2] <- base_count["C"]
    base_count_mat[n, 3] <- base_count["G"]
    base_count_mat[n, 4] <- base_count["T"]
  }
  base_count_mat[is.na(base_count_mat)] <- 0

  base_ratio_mat <- base_count_mat / rowSums(base_count_mat)
  pwm <- makePWM(t(base_ratio_mat))
 
  pdf(paste("../output/motif_diff/donor_disruption_", pos, "_normal.pdf", sep =""))
  seqLogo(pwm, ic.scale = TRUE, xaxis = FALSE, yaxis = FALSE)
  dev.off()
 
  pdf(paste("../output/motif_diff_ns/donor_disruption_", pos, "_normal.pdf", sep =""))
  seqLogo(pwm, ic.scale = FALSE, xaxis = FALSE, yaxis = FALSE)
  dev.off()
 
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "donor", Rel_Start_Motif == pos) %>%
    filter(GenomonSplicingMutation != "---")
  
  base_count_mat <- matrix(0, 9, 4)
  for(n in 1:9) {
    base_count <- table(substring(motif_sub$Motif_Seq, n + 0, n + 0))
    base_count_mat[n, 1] <- base_count["A"]
    base_count_mat[n, 2] <- base_count["C"]
    base_count_mat[n, 3] <- base_count["G"]
    base_count_mat[n, 4] <- base_count["T"]
  }
  base_count_mat[is.na(base_count_mat)] <- 0
  
  base_ratio_mat <- base_count_mat / rowSums(base_count_mat)
  pwm <- makePWM(t(base_ratio_mat))
  
  pdf(paste("../output/motif_diff/donor_disruption_", pos, "_abnormal.pdf", sep =""))
  seqLogo(pwm, ic.scale = TRUE, xaxis = FALSE, yaxis = FALSE)
  dev.off()

  pdf(paste("../output/motif_diff_ns/donor_disruption_", pos, "_abnormal.pdf", sep =""))
  seqLogo(pwm, ic.scale = FALSE, xaxis = FALSE, yaxis = FALSE)
  dev.off()

}



for(pos in 1:7) {
  
  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "acceptor", Rel_Start_Motif == pos) %>%
    filter(GenomonSplicingMutation == "---") 
  
  base_count_mat <- matrix(0, 7, 4)
  for(n in 1:7) {
    base_count <- table(substring(motif_sub$Motif_Seq, n + 0, n + 0))
    base_count_mat[n, 1] <- base_count["A"]
    base_count_mat[n, 2] <- base_count["C"]
    base_count_mat[n, 3] <- base_count["G"]
    base_count_mat[n, 4] <- base_count["T"]
  }
  base_count_mat[is.na(base_count_mat)] <- 0
  
  base_ratio_mat <- base_count_mat / rowSums(base_count_mat)
  pwm <- makePWM(t(base_ratio_mat))
  
  pdf(paste("../output/motif_diff/acceptor_disruption_", pos, "_normal.pdf", sep =""))
  seqLogo(pwm, ic.scale = TRUE, xaxis = FALSE, yaxis = FALSE)
  dev.off()

  pdf(paste("../output/motif_diff_ns/acceptor_disruption_", pos, "_normal.pdf", sep =""))
  seqLogo(pwm, ic.scale = FALSE, xaxis = FALSE, yaxis = FALSE)
  dev.off()


  motif_sub <- splicing_mutation %>% 
    filter(Type_Motif == "acceptor", Rel_Start_Motif == pos) %>%
    filter(GenomonSplicingMutation != "---") %>%
    filter(FKPM >= 1.0)
  
  base_count_mat <- matrix(0, 7, 4)
  for(n in 1:7) {
    base_count <- table(substring(motif_sub$Motif_Seq, n + 0, n + 0))
    base_count_mat[n, 1] <- base_count["A"]
    base_count_mat[n, 2] <- base_count["C"]
    base_count_mat[n, 3] <- base_count["G"]
    base_count_mat[n, 4] <- base_count["T"]
  }
  base_count_mat[is.na(base_count_mat)] <- 0
  
  base_ratio_mat <- base_count_mat / rowSums(base_count_mat)
  pwm <- makePWM(t(base_ratio_mat))
  
  pdf(paste("../output/motif_diff/acceptor_disruption_", pos, "_abnormal.pdf", sep =""))
  seqLogo(pwm, ic.scale = TRUE, xaxis = FALSE, yaxis = FALSE)
  dev.off()
  
  pdf(paste("../output/motif_diff_ns/acceptor_disruption_", pos, "_abnormal.pdf", sep =""))
  seqLogo(pwm, ic.scale = FALSE, xaxis = FALSE, yaxis = FALSE)
  dev.off()

}





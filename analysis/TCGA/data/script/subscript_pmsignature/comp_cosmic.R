library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
ctype <- args[1]
sig_num <- args[2]
input_file <- args[3]
output_file <- args[4]
cosmic_file <- args[5]

cosmic_sig_tmp <- read.table(cosmic_file, header = TRUE, sep = "\t") %>%
  arrange(Substitution.Type, Trinucleotide)

cosmic_sig <- t(cosmic_sig_tmp[1:nrow(cosmic_sig_tmp),4:35])

convertSignatureMatrixToVector <- function(Fmat, fdim) {
  
  M <- prod(fdim)
  Fvec <- rep(1, M)
  
  temp1 <- 1
  temp2 <- 1
  for (i in 1:length(fdim)) {
    temp1 <- temp1 * fdim[i]
    divInd <- (1:M - 1) %% temp1 + 1
    for (j in 1:fdim[i]) {
      targetInd <- divInd > temp2 * (j - 1) & divInd <= temp2 * j
      Fvec[targetInd] <- Fvec[targetInd] * Fmat[i,j]
    }
    temp2 <- temp2 * fdim[i]
  }
  
  return(Fvec)
}


load(input_file)


sig_ind <- c()
cosm_ind <- c()
cor_value <- c()
for (i in 1:(as.numeric(sig_num) - 1)) {
  p <- convertSignatureMatrixToVector(resultForSave[[1]]@signatureFeatureDistribution[i,c(4, 3, 1),], c(4, 4, 6))
  
  temp_cosm <- c()
  for (j in 1:32) {
    tcr <- cor(p, cosmic_sig[j,])
    if (tcr >= 0.7) {
      # temp_cosm <- c(temp_cosm, paste("cosmic_signatrue_", j, sep = ""))
      sig_ind <- c(sig_ind, i)
      cosm_ind <- c(cosm_ind, j)
      cor_value <-c (cor_value, tcr)
      
      # cosm_sig <- rbind(cosm_sig, c(sig_ind = i, coms_ind = j, cor = tcr))
    }
  }

}

cosm_sig <- data.frame(sig_ind = sig_ind, cosm_ind = cosm_ind, cor_value = cor_value)


write.table(cosm_sig, output_file, 
            row.names = FALSE, quote = FALSE, sep = "\t")



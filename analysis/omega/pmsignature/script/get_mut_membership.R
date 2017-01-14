library(pmsignature)

args <- commandArgs(trailingOnly = TRUE)
ctype <- args[1]
sig_num <- args[2]

G <- readMPFile(paste("../output/", ctype, "/pmsignature/pmsig.input.txt", sep = ""), numBases = 5, trDir = TRUE)

load(paste("../output/", ctype, "/pmsignature/ind.", sig_num, ".Rdata", sep = ""))

Q <- getMutMembership(G, resultForSave[[1]])

write.table(Q, paste("../output/", ctype, "/mapping/mut_membership.txt", sep = ""),
            row.names = FALSE, quote = FALSE, sep = "\t")


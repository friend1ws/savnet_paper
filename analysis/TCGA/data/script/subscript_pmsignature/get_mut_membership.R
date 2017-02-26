library(pmsignature)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
pmsig_result_file <- args[2]
output_file <- args[3]

G <- readMPFile(input_file, numBases = 5, trDir = TRUE)

load(pmsig_result_file)

Q <- getMutMembership(G, resultForSave[[1]])

write.table(Q, output_file,
            row.names = FALSE, quote = FALSE, sep = "\t")


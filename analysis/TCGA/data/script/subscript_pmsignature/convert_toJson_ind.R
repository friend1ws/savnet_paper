#$ R --vanilla --slave --args *.Rdata {output.json} < convert_toJson_ind.R

library(rjson)

args <- commandArgs(trailingOnly = T)
mutation_count <- list()

if (is.na(args[1])) {
    output = "output.json"
    container <- Param
} else {
    load(args[1])
    output = args[2]
    container <- resultForSave[[1]]
    if (length(resultForSave) == 3) {
        mutation_count <- resultForSave[[3]];
    }
}

cut_digits <- function(x, digits=3) {
    return (floor(x * 10^digits) /10^digits)
}
# id
ids <- list()
for (i in 1:length(container@sampleList)) {
    ids[i] <- container@sampleList[i]
}

# mutation
samples <- container@sampleSignatureDistribution
mutation <- list()
x <- 1
for (i in 1:length(samples[,1])) {
    for (j in 1:length(samples[i,])) {
        mutation[[x]] <- c(i-1, j-1, cut_digits(samples[i,j]))
        x <- x+1
    }
}

# ref, alt, strand
ref <- list()
alt <- list()
strand <- list()
for (i in 1:(container@signatureNum-1)) {
    sig <- container@signatureFeatureDistribution[i,,]
    
    r1 <- c(cut_digits(sig[2,1]),cut_digits(sig[2,2]),cut_digits(sig[2,3]),cut_digits(sig[2,4]))
    r2 <- c(cut_digits(sig[3,1]),cut_digits(sig[3,2]),cut_digits(sig[3,3]),cut_digits(sig[3,4]))
    r4 <- c(cut_digits(sig[4,1]),cut_digits(sig[4,2]),cut_digits(sig[4,3]),cut_digits(sig[4,4]))
    r5 <- c(cut_digits(sig[5,1]),cut_digits(sig[5,2]),cut_digits(sig[5,3]),cut_digits(sig[5,4]))
    r3 <- c(0, cut_digits(sig[1,1]+sig[1,2]+sig[1,3]), 0, cut_digits(sig[1,4]+sig[1,5]+sig[1,6]))
    
    ref[[i]] <- list(r1,r2,r3,r4,r5)
    
    a1 <- c(0,0,0,0)
    a2 <- c(cut_digits(sig[1,1]),0,cut_digits(sig[1,2]),cut_digits(sig[1,3]))
    a3 <- c(0,0,0,0)
    a4 <- c(cut_digits(sig[1,4]),cut_digits(sig[1,5]),cut_digits(sig[1,6]),0)
    
    alt[[i]] <- list(a1,a2,a3,a4)
    
    if (length(sig[,1]) == 6) {
        strand[[i]] <- c(cut_digits(sig[6,1]), cut_digits(sig[6,2]))
    } else {
        strand[[i]] <- list()
    }
}

# write
write(toJSON(list( id = ids, ref = ref, alt = alt, strand = strand, mutation = mutation, mutation_count = mutation_count)), output)

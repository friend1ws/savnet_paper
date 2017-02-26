#$ R --vanilla --slave --args *.Rdata {output.json} < convert_toJson_full.R

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

cut_digits <- function(x, digits=4) {
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

# signature
toList <- function(x) {
    li <- c()
    for (i in 1:length(x)) {
        #li[i] <- cut_digits(x[i]*1000000)
        li[i] <- cut_digits(x[i])
    }
    return (li)
}

sig_json <- list()
sec_num = length(container@signatureFeatureDistribution[1,,])/6
for (i in 1:(container@signatureNum-1)) {
    sig <- container@signatureFeatureDistribution[i,,]
    
    r1 <- toList(sig[1:sec_num])
    r2 <- toList(sig[(sec_num+1):(sec_num*2)])
    r3 <- toList(sig[(sec_num*2+1):(sec_num*3)])
    r4 <- toList(sig[(sec_num*3+1):(sec_num*4)])
    r5 <- toList(sig[(sec_num*4+1):(sec_num*5)])
    r6 <- toList(sig[(sec_num*5+1):(sec_num*6)])
    
    sig_json[[i]] <- list(r1,r2,r3,r4,r5,r6)
}

# write
write(toJSON(list( id = ids, signature = sig_json, mutation = mutation, mutation_count = mutation_count)), output)

#! /usr/bin/env R
#
# ${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${OUTPUTFILE} ${SIGNUM} ${TRDIRFLAG} ${TRIALNUM} ${BS_GENOME} ${TXDB_TRANSCRIPT} < run_pmsignature_ind.R
#
# ### params ###
# ${INPUTFILE}  : the path to intput mutation file. (see https://github.com/friend1ws/pmsignature)
# ${OUTPUTFILE} : the path to output file (.RData)
# ${SIGNUM}     : the number of mutation signature. (3 ~ 6)
# ${TRDIRFLAG}  : TRUE
# ${TRIALNUM}   : 10
# ${BS_GENOME}       : "BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19"
# ${TXDB_TRANSCRIPT} : "TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene"
#

library(pmsignature);

inputFile <- commandArgs()[5];
outputFile <- commandArgs()[6];
sigNum <- as.numeric(commandArgs()[7]);
trDirFlag <- as.logical(commandArgs()[8]);
trialNum <- as.numeric(commandArgs()[9]);

bs_genome <- "BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19"
txdb <- "TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene"

if (length(commandArgs()) >= 11) {
    bs_genome <- commandArgs()[10];
    txdb <- commandArgs()[11];
}

G <- readMPFile(inputFile, numBases = 5, trDir = trDirFlag, bs_genome = eval(parse(text=bs_genome)), txdb_transcript = eval(parse(text=txdb)));

BG_prob <- readBGFile(G);
Param <- getPMSignature(G, K = sigNum , BG = BG_prob, numInit = trialNum);

Boot <- bootPMSignature(G, Param0 = Param, bootNum = 100, BG = BG_prob);

mutation_count <- list()
for (i in 1:length(G@sampleList)) {
    mutation_count[i] <- length(which(G@mutationPosition==i));
}

resultForSave <- list(Param, Boot, mutation_count);
save(resultForSave, file=outputFile);


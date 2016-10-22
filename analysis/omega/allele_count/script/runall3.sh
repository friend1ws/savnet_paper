#! /bin/bash
#$ -S /bin/bash
#$ -cwd

python summarize_allele_count.py \
    ../../matome/omega.genomon_splicing_mutation.result.txt \
    ../output/ \
    ../output/omega.splicing_mutation.info.txt

Rscript splicing_mut_count.R


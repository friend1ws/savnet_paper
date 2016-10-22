#! /bin/bash
#$ -S /bin/bash

python add_allele_count.py \
    ../../matome/omega.genomon_splicing_mutation.result.txt \
    ../output/ \
    ../output/omega.genomon_splicing_mutation.allele_count.txt

Rscript IR_VAF.R


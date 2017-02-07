#! /bin/bash
#$ -S /bin/bash
#$ -cwd

python summarize_allele_count.py \
    ../../matome/omega.genomon_splicing_mutation.result.txt \
    ../output/ \
    /home/w3varann/database/GRCh37/GRCh37.fa \
    ../output/omega.splicing_mutation.info.txt

Rscript splicing_mut_count.R

Rscript mes_hb_diff_score2.R

Rscript motif_dist3.R 

Rscript motif_logo_check2.R 
 

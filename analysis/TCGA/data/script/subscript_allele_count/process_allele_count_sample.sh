#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

bam_file=$1
mut_file=$2
out_file=$3

reference=../../../db/GRCh37/GRCh37.fa

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yshira/bin/Complete-Striped-Smith-Waterman-Library/src


intron_retention_utils allele_count --grc \
    --donor_size 3,6 \
    --acceptor_size 6,1 \
     ${bam_file} ${mut_file} ${out_file} ${reference}


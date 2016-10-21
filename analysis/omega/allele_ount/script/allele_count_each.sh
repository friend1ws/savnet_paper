#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

bam_file=$1
mut_file=$2
out_file=$3

reference=/home/w3varann/database/GRCh37/GRCh37.fa
ref_gene=/home/yshira/mysoftware/intron_retention_utils/resource/refGene.txt.gz

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yshira/bin/Complete-Striped-Smith-Waterman-Library/src


intron_retention_utils allele_count \
    --chr_name_list /home/yshira/mysoftware/intron_retention_utils/resource/ucsc2grch.txt \
     ${bam_file} ${mut_file} ${out_file} ${reference} ${ref_gene}


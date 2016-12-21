#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1
SF_GENE=$2

Rscript perform_t_test.R \
    ../../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt \
    ../../gsm_out/d2.6_a8.1_8_ka/${CTYPE}/${CTYPE} \
    ../output/${CTYPE} \
    ${CTYPE} \
    ${SF_GENE}

 

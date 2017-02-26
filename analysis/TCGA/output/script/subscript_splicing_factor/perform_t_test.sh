#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1
SF_GENE=$2

Rscript subscript_splicing_factor/perform_t_test.R \
    ../../data/input_list/${CTYPE}.mut_SJ_IR_list.txt \
    ../gsm_out/d3.6_a6.1_8_ka/${CTYPE}/${CTYPE} \
    ../splicing_factor/${CTYPE} \
    ${CTYPE} \
    ${SF_GENE} \
    ../splicing_factor/sample_list_with_sfm.txt

 

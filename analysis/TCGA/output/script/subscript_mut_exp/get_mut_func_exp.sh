#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1
OUTPUT=$2

python subscript_mut_exp/get_mut_func_exp.py \
    ../../data/input_list/${CTYPE}.mut_SJ_IR_list.txt \
    ../savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../../data/expression/${CTYPE} \
    ${OUTPUT}



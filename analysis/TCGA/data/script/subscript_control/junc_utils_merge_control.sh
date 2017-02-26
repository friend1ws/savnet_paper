#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

export LANG=C

INPUT_LIST=$1
OUTPUT=$2
SAMPLE_NUM=$3

echo "junc_utils merge_control --read_num_thres 2 --keep_annotated --sample_num_thres ${SAMPLE_NUM} ${INPUT_LIST} ${OUTPUT}"
junc_utils merge_control --read_num_thres 2 --keep_annotated --sample_num_thres ${SAMPLE_NUM} ${INPUT_LIST} ${OUTPUT}


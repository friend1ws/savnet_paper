#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

export LANG=C

INPUT_LIST=$1
OUTPUT=$2
SAMPLE_NUM=$3

echo "intron_retention_utils merge_control ${INPUT_LIST} ${OUTPUT} --sample_num_thres ${SAMPLE_NUM}"
intron_retention_utils merge_control ${INPUT_LIST} ${OUTPUT} --sample_num_thres ${SAMPLE_NUM}
 

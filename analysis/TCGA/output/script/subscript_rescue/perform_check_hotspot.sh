#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

export PATH=/usr/local/package/r/current3/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/package/r/current3/lib64/R/lib:$LD_LIBRARY_PATH

export PATH=/usr/local/package/r/current3/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/package/r/current3/lib64/R/lib:$LD_LIBRARY_PATH

TUMOR_BAM=$1
NORMAL_BAM=$2
ANNOT_FILE=$3
OUTPUT=$4

python subscript_rescue/check_hotspot.py ${TUMOR_BAM} ${NORMAL_BAM} ${ANNOT_FILE} ${OUTPUT}



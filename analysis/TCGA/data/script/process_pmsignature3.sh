#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

export PATH=/usr/local/package/r/current3/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/package/r/current3/lib64/R/lib:$LD_LIBRARY_PATH
# export R_LIBS=/home/w3varann/.genomon_local/genomon_pipeline-2.5.0/R3-library
export R_LIBS=/home/yshira/.R
export R_PATH=/usr/local/package/r/current3/bin

while read line 
do

    CTYPE=`echo $line | cut -f 1 -d ' '`
    SIG_NUM=`echo $line | cut -f 2 -d ' '`

    echo "qsub subscript_pmsignature/mapping_part.sh ${CTYPE} ${SIG_NUM}"
    qsub subscript_pmsignature/mapping_part.sh ${CTYPE} ${SIG_NUM}


done < ../pmsignature/db/manual_sig_num.txt

 


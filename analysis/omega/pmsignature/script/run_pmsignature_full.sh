#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

export PATH=/usr/local/package/r/current3/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/package/r/current3/lib64/R/lib:$LD_LIBRARY_PATH
# export R_LIBS=/home/w3varann/.genomon_local/genomon_pipeline-2.5.0/R3-library
export R_LIBS=/home/yshira/.R
export R_PATH=/usr/local/package/r/current3/bin

CTYPE=$1
SIG_NUM=$2

$R_PATH/R --vanilla --slave --args ../output/${CTYPE}/pmsignature/pmsig.input.txt ../output/${CTYPE}/pmsignature/full.${SIG_NUM}.Rdata ${SIG_NUM} F 10 < run_pmsignature_full.R
if [ $? -ne 0 ]
then
    echo pmsignature terminated abnormally.
    echo "{'id':[],'signature':[],'mutation':[]}" > ../output/${CTYPE}/pmsignature/pmsignature.full.result.${SIG_NUM}.json
    exit 0
fi
$R_PATH/R --vanilla --slave --args ../output/${CTYPE}/pmsignature/full.${SIG_NUM}.Rdata ../output/${CTYPE}/pmsignature/pmsignature.full.result.${SIG_NUM}.json < convert_toJson_full.R
 

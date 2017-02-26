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

$R_PATH/R --vanilla --slave --args ../pmsignature/${CTYPE}/pmsignature/pmsig.input.txt ../pmsignature/${CTYPE}/pmsignature/ind.${SIG_NUM}.Rdata ${SIG_NUM} T 10 < subscript_pmsignature/run_pmsignature_ind.R

$R_PATH/R --vanilla --slave --args ../pmsignature/${CTYPE}/pmsignature/ind.${SIG_NUM}.Rdata ../pmsignature/${CTYPE}/pmsignature/pmsignature.ind.result.${SIG_NUM}.json < subscript_pmsignature/convert_toJson_ind.R
 

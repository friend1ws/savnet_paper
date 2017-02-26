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

mkdir -p ../pmsignature/${CTYPE}/mapping

echo "$R_PATH/R --vanilla --slave --args ../pmsignature/${CTYPE}/pmsignature/pmsig.input.txt ../pmsignature/${CTYPE}/pmsignature/ind.${SIG_NUM}.Rdata ../pmsignature/${CTYPE}/mapping/mut_membership.txt < subscript_pmsignature/get_mut_membership.R"
$R_PATH/R --vanilla --slave --args ../pmsignature/${CTYPE}/pmsignature/pmsig.input.txt ../pmsignature/${CTYPE}/pmsignature/ind.${SIG_NUM}.Rdata ../pmsignature/${CTYPE}/mapping/mut_membership.txt < subscript_pmsignature/get_mut_membership.R

echo "$R_PATH/R --vanilla --slave --args ${CTYPE} ${SIG_NUM} ../pmsignature/${CTYPE}/pmsignature/ind.${SIG_NUM}.Rdata ../pmsignature/${CTYPE}/mapping/cosm_sig.txt ../pmsignature/db/signatures_probabilities.proc.txt < subscript_pmsignature/comp_cosmic.R"
$R_PATH/R --vanilla --slave --args ${CTYPE} ${SIG_NUM} ../pmsignature/${CTYPE}/pmsignature/ind.${SIG_NUM}.Rdata ../pmsignature/${CTYPE}/mapping/cosm_sig.txt ../pmsignature/db/signatures_probabilities.proc.txt < subscript_pmsignature/comp_cosmic.R



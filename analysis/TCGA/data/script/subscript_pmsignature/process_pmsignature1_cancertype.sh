#! /bin/bash

CTYPE=$1

if [ ! -d ../pmsignature/${CTYPE}/pmsignature ]
then
    mkdir -p ../pmsignature/${CTYPE}/pmsignature
fi

echo "python subscript_pmsignature/make_pmsignature_input.py ../input_list/${CTYPE}.mut_SJ_IR_list.txt ../pmsignature/${CTYPE}/pmsignature/pmsig.input.txt"
python subscript_pmsignature/make_pmsignature_input.py ../input_list/${CTYPE}.mut_SJ_IR_list.txt ../pmsignature/${CTYPE}/pmsignature/pmsig.input.txt

for i in `seq 2 8`
do 
    echo "qsub subscript_pmsignature/run_pmsignature_ind.sh ${CTYPE} ${i}" 
    qsub subscript_pmsignature/run_pmsignature_ind.sh ${CTYPE} ${i}

    echo "qsub subscript_pmsignature/run_pmsignature_full.sh ${CTYPE} ${i}"
    qsub subscript_pmsignature/run_pmsignature_full.sh ${CTYPE} ${i}
done



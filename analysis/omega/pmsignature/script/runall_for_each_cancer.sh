#! /bin/bash

CTYPE=$1

if [ ! -d ../output/${CTYPE}/pmsignature ]
then
    mkdir -p ../output/${CTYPE}/pmsignature
fi

echo "python make_pmsignature_input.py ../../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt ../output/${CTYPE}/pmsignature/pmsig.input.txt"
python make_pmsignature_input.py ../../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt ../output/${CTYPE}/pmsignature/pmsig.input.txt

for i in `seq 2 8`
do 
    echo "qsub run_pmsignature_ind.sh ${CTYPE} ${i}" 
    qsub run_pmsignature_ind.sh ${CTYPE} ${i}

    echo "qsub run_pmsignature_full.sh ${CTYPE} ${i}"
    qsub run_pmsignature_full.sh ${CTYPE} ${i}
done


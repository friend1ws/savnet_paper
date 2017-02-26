#! /bin/bash

while read CTYPE
do

    echo "qsub subscript_savnet/execute_savnet_part.sh ${CTYPE} 3 6 6 1 8 1"
    qsub subscript_savnet/execute_savnet_part.sh ${CTYPE} 3 6 6 1 8 1

    echo "qsub subscript_savnet/execute_savnet_part.sh ${CTYPE} 5 15 15 5 8 1"
    qsub subscript_savnet/execute_savnet_part.sh ${CTYPE} 5 15 15 5 8 1

done < ../../data/input_list/cancer_type_list.txt




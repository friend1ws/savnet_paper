#! /bin/bash

while read CTYPE
do
    echo "qsub GSM_part.sh ${CTYPE}" 
    qsub GSM_part.sh ${CTYPE}
done < cancer_type_list.txt


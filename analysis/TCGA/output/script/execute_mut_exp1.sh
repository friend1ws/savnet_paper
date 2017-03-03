#! /bin/bash


if [ ! -d  ../mut_func_exp ]
then
    mkdir ../mut_func_exp
fi

while read CTYPE
do
    echo "qsub subscript_mut_exp/get_mut_func_exp.sh ${CTYPE} ../mut_func_exp/${CTYPE}.mut_func_exp.txt"
    qsub subscript_mut_exp/get_mut_func_exp.sh ${CTYPE} ../mut_func_exp/${CTYPE}.mut_func_exp.txt

done < ../../data/input_list/cancer_type_list.txt




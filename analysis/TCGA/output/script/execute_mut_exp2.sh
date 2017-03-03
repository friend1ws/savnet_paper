#! /bin/bash

echo -e "Cancer_Type\tSample_Name\tGene_Symbol\tFunc_Class\tIs_Inframe\tIs_GSM\tFPKM\tFPKM_normalized\tFPKM_mean\tFPKM_std" > ../mut_func_exp/TCGA.mut_func_exp.merged.txt

while read CTYPE
do
    echo ${CTYPE}
    while read line
    do
        echo -e "${CTYPE}\t${line}" >> ../mut_func_exp/TCGA.mut_func_exp.merged.txt
    done < ../mut_func_exp/${CTYPE}.mut_func_exp.txt
done < ../../data/input_list/cancer_type_list.txt



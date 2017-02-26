#! /bin/bash

while read CTYPE
do

    python subscript/make_savnet_input_cancertype.py ../input_list/${CTYPE}.mut_SJ_IR_list.txt ../mutation/${CTYPE} ../junction/${CTYPE} ../intron_retention/${CTYPE} ../qc/${CTYPE}

done < ../input_list/cancer_type_list.txt


#! /bin/bash

while read CTYPE 
do
    echo "bash subscript_pmsignature/paplot.sh ${CTYPE}"
    bash subscript_pmsignature/paplot.sh ${CTYPE}

done < ../input_list/cancer_type_list.txt



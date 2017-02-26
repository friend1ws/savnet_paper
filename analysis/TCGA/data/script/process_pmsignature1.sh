#! /bin/bash

while read ctype
do
    bash subscript_pmsignature/process_pmsignature1_cancertype.sh ${ctype}
done < ../input_list/cancer_type_list.txt



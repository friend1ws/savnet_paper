#! /bin/bash

while read ctype
do
    echo "bash paplot.sh ${ctype}"
    bash paplot.sh ${ctype}
    
done < cancer_type_list.txt



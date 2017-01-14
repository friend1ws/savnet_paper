#! /bin/bash

while read ctype
do
    bash runall_for_each_cancer.sh ${ctype}
done < cancer_type_list.txt



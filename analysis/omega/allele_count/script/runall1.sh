#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

while read ctype
do
    qsub runall1_for_each_cancer.sh ${ctype}
    sleep 20
    # sleep 0
done < cancer_type_list.txt
 

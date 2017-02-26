#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

while read ctype
do
    bash subscript_allele_count/process_allele_count_cancertype.sh ${ctype}
    sleep 20
    # sleep 0
done < ../input_list/cancer_type_list.txt 

 

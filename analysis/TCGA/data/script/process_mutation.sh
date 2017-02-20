#! /bin/bash
#$ -S /bin/bash
#$ -cwd

if [ ! -d log ]
then
    mkdir log
fi

annot_utils boundary --donor_size 0,2 --acceptor_size 2,0 --grc ../mutation/canonical.bed.gz

annot_utils boundary --donor_size 5,15 --acceptor_size 15,5 --grc ../mutation/noncanonical.bed.gz



while read ctype
do
    echo "qsub remove_dup_cancer_type.sh ${ctype}"
    qsub remove_dup_cancer_type.sh ${ctype}
done < ../input_list/cancer_type_list.txt



#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1

if [ ! -d ../allele_count/${CTYPE} ]
then
    mkdir -p ../allele_count/${CTYPE}
fi

python subscript_allele_count/generate_mut_bam_list.py \
    ../allele_count/${CTYPE}/mut_bam_list.txt \
    ../mutation/${CTYPE} \
    /home/omega3/omega_rna/star/${CTYPE} \
    /home/eva/rawdata/tcga_rna_prev/single/output/${CTYPE}

while read line
do
    sample=`echo ${line} | cut -d ' ' -f 1`
    mut_file=`echo ${line} | cut -d ' ' -f 2`
    bam_file=`echo ${line} | cut -d ' ' -f 3`

    # if [ ! -s ../allele_count/${CTYPE}/${sample}.allele_count.txt ]
    # then
        echo "qsub subscript_allele_count/process_allele_count_sample.sh ${bam_file} ${mut_file} ../allele_count/${CTYPE}/${sample}.allele_count.txt"
        qsub subscript_allele_count/process_allele_count_sample.sh ${bam_file} ${mut_file} ../allele_count/${CTYPE}/${sample}.allele_count.txt
    # fi

done < ../allele_count/${CTYPE}/mut_bam_list.txt




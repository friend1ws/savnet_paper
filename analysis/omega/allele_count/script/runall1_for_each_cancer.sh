#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1

if [ ! -d ../output/${CTYPE} ]
then
    mkdir -p ../output/${CTYPE}
fi

python generate_mut_bam_list.py \
    ../output/${CTYPE}/mut_bam_list.txt \
    /home/yshira/mypaper/mutrans_paper/analysis/omega/mutation/output/${CTYPE} \
    /home/omega3/omega_rna/star/${CTYPE}

while read line
do
    sample=`echo ${line} | cut -d ' ' -f 1`
    mut_file=`echo ${line} | cut -d ' ' -f 2`
    bam_file=`echo ${line} | cut -d ' ' -f 3`

    # if [ ! -s ../output/${CTYPE}/${sample}.allele_count.txt ]
    # then
        echo "qsub allele_count_each.sh ${bam_file} ${mut_file} ../output/${CTYPE}/${sample}.allele_count.txt"
        qsub allele_count_each.sh ${bam_file} ${mut_file} ../output/${CTYPE}/${sample}.allele_count.txt
        # echo "bash allele_count_each.sh ${bam_file} ${mut_file} ../output/${CTYPE}/${sample}.allele_count.txt"
        # bash allele_count_each.sh ${bam_file} ${mut_file} ../output/${CTYPE}/${sample}.allele_count.txt
    # fi

done < ../output/${CTYPE}/mut_bam_list.txt



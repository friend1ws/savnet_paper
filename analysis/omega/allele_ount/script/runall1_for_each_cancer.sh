#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1

if [ ! -d ../output/${CTYPE} ]
then
    mkdir -p ../output/${CTYPE}
fi

python generate_omega_mut_list.py \
    ../output/${CTYPE}/mut_bam_list.txt \
    /home/kchiba/work_directory/work_hotspot/black_list_output_min10/${CTYPE} \
    /home/omega3/omega_rna/star/${CTYPE}

while read line
do
    sample=`echo ${line} | cut -d ' ' -f 1`
    mut_file=`echo ${line} | cut -d ' ' -f 2`
    bam_file=`echo ${line} | cut -d ' ' -f 3`

    echo "qsub allele_count_each.sh ${bam_file} ${mut_file} ../output/${CTYPE}/${sample}.allele_count.txt"
    qsub allele_count_each.sh ${bam_file} ${mut_file} ../output/${CTYPE}/${sample}.allele_count.txt

done < ../output/${CTYPE}/mut_bam_list.txt



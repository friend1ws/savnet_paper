#! /bin/bash
#$ -S /bin/bash
#$ -cwd

if [ ! -d ../rescue ]
then
    mkdir -p ../rescue
fi

# First, list up the SAVs detected by SAVNet
python subscript_rescue/generate_vcf_list.py ../savnet_out/d3.6_a6.1_8_ka/TCGA.genomon_splicing_mutation.result.txt | sort -k1,1 -k2,2n > ../rescue/sav.anno.txt


# Then for each sample, extract the somatic variants that satisfy the following criteria
# 1. Number of variant supporting read is not less than 3
# 2. variant allele frequency (VAF)of tumor is not less than 0.1
# 3. the ratio between VAF of tumor and control is not less than 10

while read CTYPE 
do

    if [ ! -d ../rescue/${CTYPE} ]
    then
        mkdir -p ../rescue/${CTYPE}
    fi

    for ID_T in `cut -f 1 ../../data/input_list/${CTYPE}.mut_SJ_IR_list.txt | tail -n +2`
    do
        ID_N=${ID_T%%01}
        ID_N=${ID_N%%06}1

        TUMOR_BAM=/home/omega3/omega_project/genomon/${CTYPE}/bam/${ID_T}/${ID_T}.markdup.bam
        NORMAL_BAM=`ls /home/omega3/omega_project/genomon/${CTYPE}/bam/${ID_N}*/${ID_N}*.markdup.bam`

        if [ ! -s ../rescue/${CTYPE}/${ID_T}.check_hotspot.txt ]
        then
            echo "qsub subscript_rescue/perform_check_hotspot.sh ${TUMOR_BAM} ${NORMAL_BAM} ../rescue/sav.anno.txt ../rescue/${CTYPE}/${ID_T}.check_hotspot.txt"
            qsub subscript_rescue/perform_check_hotspot.sh ${TUMOR_BAM} ${NORMAL_BAM} ../rescue/sav.anno.txt ../rescue/${CTYPE}/${ID_T}.check_hotspot.txt
        fi

    done

    sleep 900

done < ../../data/input_list/cancer_type_list.txt

# done < ACC.txt

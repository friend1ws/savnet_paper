#! /bin/bash

if [ ! -d ../gsm_out/gsm_file_list ]
then
    mkdir -p ../gsm_out/gsm_file_list
fi

while read CTYPE
do
    echo $CTYPE

    python generate_omega_mut_list.py \
        ../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt \
        /home/yshira/mypaper/mutrans_paper/analysis/omega/mutation/output/${CTYPE} \
        /home/eva/genomon_out/rna_2_4_0/TCGA/${CTYPE}/star \
        /home/eva/rawdata/tcga_rna_prev/single/output/${CTYPE}/star \
        /home/eva/genomon_out/rna_2_4_0/TCGA/${CTYPE}/intron_retention \
        /home/yshira/mypaper/mutrans_paper/analysis/omega/intron_retention_single/output/${CTYPE}

done < cancer_type_list.txt


while read CTYPE
do

    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 4 0
    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 8 0
    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 38 0   
    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 75 0
    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 1000 0

    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 4 0
    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 8 0
    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 38 0
    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 75 0
    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 1000 0

    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 4 1
    qsub GSM_part.sh ${CTYPE} 3 6 6 1 8 1
    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 38 1
    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 75 1
    # qsub GSM_part.sh ${CTYPE} 3 6 6 1 1000 1

    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 4 1
    qsub GSM_part.sh ${CTYPE} 5 15 15 5 8 1
    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 38 1
    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 75 1
    # qsub GSM_part.sh ${CTYPE} 5 15 15 5 1000 1


done < cancer_type_list.txt



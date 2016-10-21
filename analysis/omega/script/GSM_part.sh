#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1

if [ ! -d ../output/${CTYPE} ]
then
    mkdir -p ../output/${CTYPE}
fi

<<_COM_

python generate_omega_mut_list.py \
    ../output/${CTYPE}/mut_SJ_IR_list.txt \
    /home/kchiba/work_directory/work_hotspot/black_list_output_min10/${CTYPE} \
    /home/eva/genomon_out/rna_2_4_0/TCGA/${CTYPE}/star \
    /home/eva/rawdata/tcga_rna_prev/single/output/${CTYPE}/star \
    /home/eva/genomon_out/rna_2_4_0/TCGA/${CTYPE}/intron_retention \
    /home/eva/rawdata/tcga_rna_prev/single/output/${CTYPE}/intron_retention \

_COM_

genomon_splicing_mutation ../output/${CTYPE}/mut_SJ_IR_list.txt ../output/${CTYPE}/${CTYPE} /home/yshira/mysoftware/junc_utils/resource /home/w3varann/database/GRCh37/GRCh37.fa --SJ_pooled_control_file /home/yshira/project/inframe_junc/output/control.bed.gz --IR_pooled_control_file /home/yshira/project/GIR/simple_count/TCGA/output/CTRL/control.bed.gz 
 

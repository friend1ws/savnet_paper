#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1
D_E=$2
D_I=$3
A_I=$4
A_E=$5
C=$6
KA=$7

if [ $KA = 0 ]
then
    out_dir=../gsm_out/d${D_E}.${D_I}_a${A_I}.${A_E}_${C}/${CTYPE}
else
    out_dir=../gsm_out/d${D_E}.${D_I}_a${A_I}.${A_E}_${C}_ka/${CTYPE}
fi

if [ ! ${out_dir} ]
then
    mkdir -p ${out_dir}
fi

if [ $KA = 0 ]
then
    echo "genomon_splicing_mutation --grc --donor_size ${D_E},${D_I} --acceptor_size ${A_I},${A_E} ../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt ${out_dir}/${CTYPE} /home/w3varann/database/GRCh37/GRCh37.fa --SJ_pooled_control_file ../control/SJ/output/control_2_${C}.bed.gz --IR_pooled_control_file ../control/IR/output/control_${C}.bed.gz"
    genomon_splicing_mutation --grc --donor_size ${D_E},${D_I} --acceptor_size ${A_I},${A_E} ../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt ${out_dir}/${CTYPE} /home/w3varann/database/GRCh37/GRCh37.fa --SJ_pooled_control_file ../control/SJ/output/control_2_${C}.bed.gz --IR_pooled_control_file ../control/IR/output/control_${C}.bed.gz 
else
    echo "genomon_splicing_mutation --grc --keep_annotated --donor_size ${D_E},${D_I} --acceptor_size ${A_I},${A_E} ../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt ${out_dir}/${CTYPE} /home/w3varann/database/GRCh37/GRCh37.fa --SJ_pooled_control_file ../control/SJ/output/control_2_${C}.bed.gz --IR_pooled_control_file ../control/IR/output/control_${C}.bed.gz "
    genomon_splicing_mutation --grc --keep_annotated --donor_size ${D_E},${D_I} --acceptor_size ${A_I},${A_E} ../gsm_out/gsm_file_list/${CTYPE}.mut_SJ_IR_list.txt ${out_dir}/${CTYPE} /home/w3varann/database/GRCh37/GRCh37.fa --SJ_pooled_control_file ../control/SJ/output/control_2_${C}.bed.gz --IR_pooled_control_file ../control/IR/output/control_${C}.bed.gz 
fi
 

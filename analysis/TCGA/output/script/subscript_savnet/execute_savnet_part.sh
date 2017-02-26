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


out_dir=../savnet_out/d${D_E}.${D_I}_a${A_I}.${A_E}_${C}_ka/${CTYPE}

if [ ! ${out_dir} ]
then
    mkdir -p ${out_dir}
fi


echo "savnet --grc --keep_annotated --donor_size ${D_E},${D_I} --acceptor_size ${A_I},${A_E} ../../data/input_list/${CTYPE}.mut_SJ_IR_list.txt ${out_dir}/${CTYPE} ../../../db/GRCh37/GRCh37.fa --SJ_pooled_control_file ../../data/control/SJ/control_2_${C}.bed.gz --IR_pooled_control_file ../../data/control/IR/control_${C}.bed.gz "
savnet --grc --keep_annotated --donor_size ${D_E},${D_I} --acceptor_size ${A_I},${A_E} ../../data/input_list/${CTYPE}.mut_SJ_IR_list.txt ${out_dir}/${CTYPE} ../../../db/GRCh37/GRCh37.fa --SJ_pooled_control_file ../../data/control/SJ/control_2_${C}.bed.gz --IR_pooled_control_file ../../data/control/IR/control_${C}.bed.gz 
 

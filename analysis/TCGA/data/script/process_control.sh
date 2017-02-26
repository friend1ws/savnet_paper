#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

export LANG=C

if [ ! -d ../control/SJ ]
then
    mkdir -p ../control/SJ
fi

if [ ! -d ../control/IR ]
then
    mkdir -p ../control/IR
fi


echo "python python subscript_control/make_control_list_SJ.py > ../control/SJ/control_rna_list.txt"
python subscript_control/make_control_list_SJ.py > ../control/SJ/control_rna_list.txt

echo "qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_4.bed.gz 4"
qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_4.bed.gz 4

echo "qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_8.bed.gz 8"
qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_8.bed.gz 8

echo "qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_38.bed.gz 38"
qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_38.bed.gz 38

echo "qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_75.bed.gz 75"
qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_75.bed.gz 75

echo "qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_1000.bed.gz 1000"
qsub subscript_control/junc_utils_merge_control.sh ../control/SJ/control_rna_list.txt ../control/SJ/control_2_1000.bed.gz 1000


echo "python subscript_control/make_control_list_IR.py > ../control/IR/control_rna_list.txt"
python subscript_control/make_control_list_IR.py > ../control/IR/control_rna_list.txt 

echo "qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_4.bed.gz 4"
qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_4.bed.gz 4

echo "qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_8.bed.gz 8"
qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_8.bed.gz 8 

echo "qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_38.bed.gz 38"
qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_38.bed.gz 38

echo "qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_75.bed.gz 75"
qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_75.bed.gz 75

echo "qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_1000.bed.gz 1000"
qsub subscript_control/intron_retention_utils_merge_control.sh ../control/IR/control_rna_list.txt ../control/IR/control_1000.bed.gz 1000

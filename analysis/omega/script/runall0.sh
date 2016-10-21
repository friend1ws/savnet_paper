#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

if [ ! -d ../output/control ]
then
    mkdir -p ../output/control
fi

python get_control_list_SJ.py ../output/control/control_junc_list.txt

junc_utils merge_control --read_num_thres 2 --sample_num_thres 2 ../output/control/control_junc_list.txt ../output/control/control_junc.bed.gz


# python get_control_list_intron.py ../output/control/control_intron_list.txt

# intron_retention_utils merge_control ../output/control/control_intron_list.txt ../output/control/control_intron.bed.gz --sample_num_thres 5


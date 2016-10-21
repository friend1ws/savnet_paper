#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

CTYPE=$1

if [ ! -d ../output_sv/${CTYPE} ]
then
    mkdir -p ../output_sv/${CTYPE}
fi

python generate_omega_sv_list.py \
    ../output_sv/${CTYPE}/sv_SJ_IR_Chimera_list.txt \
    /home/yshira/mypaper/genomonSV_paper/analysis/omega/matome/${CTYPE} \
    /home/omega3/omega_rna/star/${CTYPE} \
    /home/yshira/project/GIR/simple_count/TCGA/output/${CTYPE} \
    /home/yshira/project/chimera/output/${CTYPE}


genomon_splicing_mutation ../output_sv/${CTYPE}/sv_SJ_IR_Chimera_list.txt ../output_sv/${CTYPE}/${CTYPE} /home/yshira/mysoftware/junc_utils/resource  /home/w3varann/database/GRCh37/GRCh37.fa --sv --SJ_pooled_control_file /home/yshira/project/inframe_junc/output/control.bed.gz --IR_pooled_control_file /home/yshira/project/GIR/simple_count/TCGA/output/CTRL/control.bed.gz --chimera_pooled_control_file /home/yshira/project/chimera/output/control.merge.bed.gz 
 

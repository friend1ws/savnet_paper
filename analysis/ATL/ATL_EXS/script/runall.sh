#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

python generate_omega_mut_list.py \
    sample_list.txt \
    /home/yshira/mypaper/genomonSV_paper/analysis/ATL/result/genomon_2_3_0/EXS/mutation \
    /home/yshira/project/ATL_RNAseq_Geno2_3_0/star \
    /home/yshira/project/GIR/simple_count/ATL/output

<<_C_

ls /home/yshira/mypaper/mutrans_paper/analysis/ATL/output/genomon_2_3_0/star/Career*/Career*.SJ.out.tab > SJ_control_list.txt
ls /home/yshira/mypaper/mutrans_paper/analysis/ATL/output/genomon_2_3_0/star/PB_CD4_*/PB_CD4_*.SJ.out.tab >> SJ_control_list.txt

junc_utils merge_control SJ_control_list.txt ../output/SJ_control.txt.gz
  
ls /home/yshira/project/GIR/simple_count/ATL/output/Career*.genomon_intron_retention.result.txt > IR_control_list.txt
ls /home/yshira/project/GIR/simple_count/ATL/output/PB_CD4_*.genomon_intron_retention.result.txt >> IR_control_list.txt

genomon_intron_retention merge_control IR_control_list.txt ../output/IR_contorl.txt.gz

_C_

genomon_splicing_mutation \
     sample_list.txt \
    ../output/ATL_EXS \
    /home/yshira/mysoftware/junc_utils/resource \
    /home/w3varann/database/GRCh37/GRCh37.fa \
    --SJ_pooled_control_file /home/yshira/project/inframe_junc/output/control.bed.gz \
     --IR_pooled_control_file /home/yshira/project/GIR/simple_count/TCGA/output/CTRL/control.bed.gz
 

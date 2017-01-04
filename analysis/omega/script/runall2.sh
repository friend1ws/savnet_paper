#! /bin/bash

if [ ! -d ../matome ]
then
    mkdir -p ../matome
fi


python subscript_matome/matome_omega.py ../gsm_out/d3.6_a6.1_8_ka ../matome/omega.genomon_splicing_mutation.result.txt /home/yshira/mysoftware/sv_utils/cancer_gene/cancer_gene.txt ../HGMD/HGMD.CS.bed.gz 3 0.05 

python subscript_matome/summarize_snv_info.py  ../matome/omega.genomon_splicing_mutation.result.txt ../matome/omega.motif_summary.txt /home/w3varann/database/GRCh37/GRCh37.fa ../spidex/hg19_spidex.bed.gz

Rscript subscript_matome/add_mes_score.R ../matome/omega.motif_summary.txt ../matome/omega.motif_summary.mes.txt

python subscript_matome/check_inframe.py ../matome/omega.genomon_splicing_mutation.result.txt ../matome/gene.inframe_summary.txt

python subscript_matome/summarize_mut_count.py ../matome/omega.mut_count.txt

python subscript_matome/check_position_fdr.py ../gsm_out/d5.15_a15.5_8_ka/ 3.0 > ../matome/position_fdr.txt

python subscript_matome/check_alt_junc.py ../matome/omega.genomon_splicing_mutation.result.txt ../matome/omega.alt_junc.txt

python subscript_matome/merge_mut.py ../matome/omega.genomon_splicing_mutation.result.txt ../gsm_out/gsm_file_list/ > ../matome/omega.mutation.merged.txt


Rscript subscript_matome/category2.R 

# Rscript subscript_matome/motif_dist.R
Rscript subscript_matome/motif_dist_creation.R

Rscript subscript_matome/multiple_events.R

Rscript subscript_matome/count_summary.R

Rscript subscript_matome/fdr_summary.R

Rscript subscript_matome/cancer_gene_summary.R 

Rscript subscript_matome/cancer_gene_ratio.R

Rscript subscript_matome/mes_hbond_spidex.fig.R

Rscript subscript_matome/alt_junc_pos.R

Rscript subscript_matome/JungEtAl_venn.R


<<_COM_

python summarize_motif_pos.py omega.genomon_splicing_mutation.result.perm.txt /home/w3varann/database/GRCh37/GRCh37.fa perm.motif_summary.txt


python check_inframe.py omega.genomon_splicing_mutation.result.txt inframe_summary.txt

_COM_


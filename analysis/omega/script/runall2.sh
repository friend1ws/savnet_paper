#! /bin/bash

if [ ! -d ../matome ]
then
    mkdir -p ../matome
fi


python subscript_matome/matome_omega.py ../output ../matome/omega.genomon_splicing_mutation.result.txt /home/yshira/mysoftware/sv_utils/cancer_gene/cancer_gene.txt 3 0.10 

python subscript_matome/summarize_motif_pos.py ../matome/omega.genomon_splicing_mutation.result.txt ../matome/omega.motif_summary.txt /home/w3varann/database/GRCh37/GRCh37.fa

python subscript_matome/check_inframe.py ../matome/omega.genomon_splicing_mutation.result.txt ../matome/gene.inframe_summary.txt

python subscript_matome/summarize_mut_count.py ../matome/omega.mut_count.txt

Rscript subscript_matome/category2.R 

Rscript subscript_matome/motif_dist.R

Rscript subscript_matome/multiple_events.R

Rscript subscript_matome/count_summary.R

Rscript subscript_matome/cancer_gene_summary.R 

<<_COM_

python summarize_motif_pos.py omega.genomon_splicing_mutation.result.perm.txt /home/w3varann/database/GRCh37/GRCh37.fa perm.motif_summary.txt


python check_inframe.py omega.genomon_splicing_mutation.result.txt inframe_summary.txt

_COM_


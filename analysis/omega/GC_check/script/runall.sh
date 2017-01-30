#! /bin/bash

python summarize_intron_info.py ../../allele_count/output/omega.splicing_mutation.info.txt ../output/omega.gsm.gc.txt /home/yshira/common/ref/GRCh37-lite_PCAWG_bwa-0.7.12/GRCh37-lite_PCAWG.fa

Rscript GC_contents_diff.R


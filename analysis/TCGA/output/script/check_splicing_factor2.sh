#! /bin/bash

python subscript_splicing_factor/check_savnet.py ../savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt > ../splicing_factor/TCGA.savnet.sf_mut.summary.txt

python subscript_splicing_factor/summarize_splicing_factor.py > ../splicing_factor/TCGA.SJ_IR.sig_summary.txt

 

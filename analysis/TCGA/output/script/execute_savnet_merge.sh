#! /bin/bash

python subscript_savnet/execute_savnet_merge.py ../savnet_out/d3.6_a6.1_8_ka ../savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt ../../../db/cancer_gene/cancer_gene.txt ../../../db/HGMD/HGMD.CS.bed.gz 3 0.05

cat ../savnet_out/d3.6_a6.1_8_ka/TCGA.savnet.result.txt \
    ../rescue/TCGA.savnet.rescued.result.txt > \
    ../rescue/TCGA.savnet.with_rescued.result.txt



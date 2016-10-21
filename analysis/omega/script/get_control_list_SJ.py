#! /usr/bin/env python

import sys, os, glob, re

output_file = sys.argv[1]

all_files = glob.glob("/home/eva/genomon_out/rna_2_4_0/TCGA/*/star/*/*.SJ.out.tab") + \
            glob.glob("/home/eva/rawdata/tcga_rna_prev/single/output/*/star/*/*.SJ.out.tab")

hout = open(output_file, 'w')
for junc_file in sorted(all_files):

    bfile = os.path.basename(junc_file)
    sample = bfile.replace(".SJ.out.tab", "")
    sample_type = re.sub(r'^0', '', sample[13:15])
    if int(sample_type) >= 10:
        print >> hout, junc_file



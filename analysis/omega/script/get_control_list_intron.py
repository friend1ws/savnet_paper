#! /usr/bin/env python

import sys, os, glob, re

output_file = sys.argv[1]

all_files = glob.glob("/home/eva/genomon_out/rna_2_4_0/TCGA/*/intron_retention/*/*.ir_simple_count.txt") + \
            glob.glob("/home/eva/rawdata/tcga_rna_prev/single/output/*/intron_retention/*/*.ir_simple_count.txt")

hout = open(output_file, 'w')
for intron_file in sorted(all_files):

    bfile = os.path.basename(intron_file)
    sample = bfile.replace(".ir_simple_count.txt", "")
    sample_type = re.sub(r'^0', '', sample[13:15])
    if int(sample_type) >= 10:
        print >> hout, intron_file



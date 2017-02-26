#! /usr/bin/env python

import glob, re, os
# junc_path = "/home/omega3/omega_rna/star"

junc_file_list = glob.glob("../junction/*/TCGA-*.SJ.out.tab")

for junc_file in sorted(junc_file_list):

    bfile = os.path.basename(junc_file)
    sample = bfile.replace(".SJ.out.tab", "")
    sample_type = re.sub(r'^0', '', sample[13:15])
    if int(sample_type) >= 10:
        # print >> hout, junc_file
        print junc_file

# hout.close()


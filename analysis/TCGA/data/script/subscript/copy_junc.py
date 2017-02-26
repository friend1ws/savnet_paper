#! /usr/bin/env python

import sys, os, glob, shutil

junc_files = glob.glob("/home/eva/genomon_out/rna_2_4_0/TCGA/*/star/*/*.SJ.out.tab") + \
                glob.glob("/home/eva/rawdata/tcga_rna_prev/single/output/*/star/*/*.SJ.out.tab")


for junc_file in sorted(junc_files):
    cancer_type = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(junc_file))))
    out_file = "../junction/" + cancer_type + "/" + os.path.basename(junc_file)

    if cancer_type in ["LAML", "OV"]: continue
    if not os.path.exists("../junction/" + cancer_type): os.makedirs("../junction/" + cancer_type)
    print cancer_type + '\t' + os.path.basename(junc_file)

    shutil.copy2(junc_file, out_file)


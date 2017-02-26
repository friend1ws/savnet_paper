#! /usr/bin/env python

import sys, os, glob, shutil

qc_files = glob.glob("/home/eva/genomon_out/rna_2_4_0/TCGA/*/star/*/*.Log.final.out") + \
                glob.glob("/home/eva/rawdata/tcga_rna_prev/single/output/*/star/*/*.Log.final.out")



for qc_file in sorted(qc_files):
    cancer_type = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(qc_file))))
    out_file = "../qc/" + cancer_type + "/" + os.path.basename(qc_file)

    if cancer_type in ["LAML", "OV"]: continue
    if not os.path.exists("../qc/" + cancer_type): os.makedirs("../qc/" + cancer_type)
    print cancer_type + '\t' + os.path.basename(qc_file)

    shutil.copy2(qc_file, out_file)


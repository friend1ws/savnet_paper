#! /usr/bin/env python

import sys, os, glob, shutil

exp_files = glob.glob("/home/eva/genomon_out/rna_2_4_0/TCGA/*/expression/*/*.sym2fkpm.txt")


for exp_file in sorted(exp_files):
    cancer_type = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(exp_file))))
    out_file = "../expression/" + cancer_type + "/" + os.path.basename(exp_file)

    if cancer_type in ["LAML", "OV"]: continue
    if not os.path.exists("../expression/" + cancer_type): os.makedirs("../expression/" + cancer_type)
    print cancer_type + '\t' + os.path.basename(exp_file)

    shutil.copy2(exp_file, out_file)


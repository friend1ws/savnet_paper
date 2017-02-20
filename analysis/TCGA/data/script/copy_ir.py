#! /usr/bin/env python

import sys, os, glob, shutil

# junc_files = glob.glob("/home/eva/genomon_out/rna_2_4_0/TCGA/*/star/*/*.SJ.out.tab") + \
#              glob.glob("/home/eva/rawdata/tcga_rna_prev/single/output/*/star/*/*.SJ.out.tab")

ir_files = glob.glob("/home/eva/genomon_out/rna_2_4_0/TCGA/*/intron_retention/*/*.ir_simple_count.txt") 

for ir_file in sorted(ir_files):
    cancer_type = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(ir_file))))
    out_file = "../intron_retention/" + cancer_type + "/" + os.path.basename(ir_file)

    if cancer_type in ["LAML", "OV"]: continue
    if not os.path.exists("../intron_retention/" + cancer_type): os.makedirs("../intron_retention/" + cancer_type)
    print cancer_type + '\t' + os.path.basename(ir_file)

    shutil.copy2(ir_file, out_file)



ir_files = glob.glob("../../../../omega/intron_retention_single/output/*/*/*.ir_simple_count.txt")

for ir_file in sorted(ir_files):
    cancer_type = os.path.basename(os.path.dirname(os.path.dirname(ir_file)))
    out_file = "../intron_retention/" + cancer_type + "/" + os.path.basename(ir_file)

    if cancer_type in ["LAML", "OV"]: continue
    if not os.path.exists("../intron_retention/" + cancer_type): os.makedirs("../intron_retention/" + cancer_type)
    print cancer_type + '\t' + os.path.basename(ir_file)
    
    shutil.copy2(ir_file, out_file)



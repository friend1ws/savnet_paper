#! /usr/bin/env python

import sys, glob, os, re

gsm_dir = sys.argv[1]
score_thres = sys.argv[2]

all_gsm_file = glob.glob(gsm_dir + "*/*.savnet.result.txt")


cancer_type2max_fdr = {}
for gsm_file in sorted(all_gsm_file):

    cancer_type = os.path.basename(os.path.dirname(gsm_file))
    with open(gsm_file, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            if float(F[header2ind["Score"]]) < float(score_thres): continue

            if cancer_type not in cancer_type2max_fdr or float(F[header2ind["Q_Value"]]) > cancer_type2max_fdr[cancer_type]:
                cancer_type2max_fdr[cancer_type] = float(F[header2ind["Q_Value"]])

print "Cancer_Type" + '\t' + "FDR"
for cancer_type in sorted(cancer_type2max_fdr):
    print cancer_type + '\t' + str(cancer_type2max_fdr[cancer_type])
        



#! /usr/bin/env python

import glob, sys, os

input_dir = sys.argv[1]
output_file = sys.argv[2]

allfiles = glob.glob(input_dir + "/*/*.genomon_splicing_mutation.result.txt")

header_order = ["Sample_Name", "SV_Key", "SV_Type", "Gene_Sybmol", "Splicing_Key", "Splicing_Class", "Is_Inframe"]

header_ind = 0
header2ind = {}
hout = open(output_file, 'w')

print >> hout, "Sample_Name" + '\t' + "Cancer_Type" + '\t' + '\t'.join(header_order[1:])

for gsm_file in sorted(allfiles):
    cancer_type = os.path.basename(os.path.dirname(gsm_file))
    with open(gsm_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i


        for line in hin:
            F = line.rstrip('\n').split('\t')
            if float(F[header2ind["Score"]]) < 3.0: continue
            print >> hout, F[header2ind["Sample_Name"]] + '\t' + cancer_type + '\t' + '\t'.join([F[header2ind[x]] for x in header_order[1:]])

hout.close()




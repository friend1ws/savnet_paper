#! /usr/bin/env python

import glob, sys, os

all_files = glob.glob("../data/gdac.broadinstitute.org_*-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt")

input_file = sys.argv[1]
output_file = sys.argv[2]

target_key = {}
with open(input_file, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = F[header2ind["Gene_Symbol"]] + '\t' + F[header2ind["Sample_Name"]]
        target_key[key] = 1


key2cna = {}
for cna_file in sorted(all_files):
    cancer_type = os.path.basename(os.path.dirname(cna_file)).replace("gdac.broadinstitute.org_", '').replace("-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0", '')
    print >> sys.stderr, "processing: " + cancer_type
    with open(cna_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        ind2sample_name = {}
        for i in range(3, len(header)):
            ind2sample_name[i] = header[i][:15]
            # print str(i) + '\t' + header[i][:15]

        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene_symbol = F[0]
            for i in range(3, len(F)):
                sample_name = ind2sample_name[i]
                key = gene_symbol + '\t' + sample_name
                if key in target_key:
                    key2cna[key] = str(F[i])
                    # print cancer_type + '\t' + key + '\t' + str(F[i])

hout = open(output_file, 'w')
with open(input_file, 'r') as hin:
    header_line = hin.readline().rstrip('\n')
    print >> hout, header_line + '\t' + "CNA"
 
    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = F[header2ind["Gene_Symbol"]] + '\t' + F[header2ind["Sample_Name"]]
        cna = key2cna[key] if key in key2cna else "NA"
        print >> hout, '\t'.join(F) + '\t' + cna

hout.close()


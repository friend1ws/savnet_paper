#! /usr/bin/env python

import sys, os, glob

annot_file = sys.argv[1]
hotspot_dir = sys.argv[2]
output_file = sys.argv[3]

mut_key2gene = {}
with open(annot_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut_key2gene['\t'.join(F[:5])] = F[5]


hotspot_files = glob.glob(hotspot_dir + "/*/*.check_hotspot.txt")
header_flag = False
hout = open(output_file, 'w')
for hfile in sorted(hotspot_files):
    sample_name = os.path.basename(hfile).replace(".check_hotspot.txt", "")
    cancer_type = os.path.basename(os.path.dirname(hfile))

    with open(hfile, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        header2ind = {}
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        if header_flag == False:
            print >> hout, "Cancer_Type" + '\t' + "Sample_Name" + '\t' + "Gene_Symbol" + '\t' + \
                    '\t'.join(header)
            header_flag = True

            
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mut_key = '\t'.join([F[header2ind[x]] for x in ["Chr", "Start", "End", "Ref", "Alt"]])
  
            if F[header2ind["Ref"]] == "-" or F[header2ind["Alt"]] == "-": continue 
            gene = mut_key2gene[mut_key]
            print >> hout, cancer_type + '\t' + sample_name + '\t' + gene + '\t' + '\t'.join(F)

 


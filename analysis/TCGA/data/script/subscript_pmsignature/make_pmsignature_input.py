#! /usr/bin/env python

import sys

gsm_input = sys.argv[1]
output_file = sys.argv[2]

hout = open(output_file, 'w')
with open(gsm_input, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        sample_name = F[header2ind["Sample_Name"]]
        mut_file = F[header2ind["Mutation_File"]]
        with open(mut_file, 'r') as hin2:
            header2ind2 = {}
            header2 = hin2.readline().rstrip('\n').split('\t')
            for (i, cname) in enumerate(header2):
                header2ind2[cname] = i

            for line2 in hin2:
                F2 = line2.rstrip('\n').split('\t')
                if F2[header2ind2["Ref"]] == "-" or F2[header2ind2["Alt"]] == "-": continue
                if F2[header2ind2["Chr"]].startswith("GL0") or F2[header2ind2["Chr"]] == "MT": continue
   
                F2[header2ind2["Chr"]] = "chr" + F2[header2ind2["Chr"]]
 
                print >> hout, sample_name + '\t' + '\t'.join([F2[header2ind2[x]] for x in ["Chr", "Start", "Ref", "Alt"]])
 



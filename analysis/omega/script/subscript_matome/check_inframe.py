#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]


key2inframe = {}
header2ind = {}

key_header = ["Cancer_Type", "Sample_Name", "Gene_Symbol", "Mutation_Key",
              "Motif_Pos", "Mutation_Type", "Is_Canonical", "LawrenceEtAl_2014",
              "Cancer_Gene_Census", "Is_Cancer_Gene"]

hout = open(output_file, 'w')
print >> hout, '\t'.join(key_header) + '\t' + "Is_Inframe"
with open(input_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    for i in range(len(header)):
        header2ind[header[i]] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = '\t'.join([F[header2ind[x]] for x in key_header])

        if key not in key2inframe: key2inframe[key] = "No"

        if F[header2ind["Is_Inframe"]] == "in-frame": key2inframe[key] = "Yes"

for key in sorted(key2inframe):
    print >> hout,key + '\t' + key2inframe[key]
                  
hout.close()                                   


#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]


key2count = {}
header2ind = {}

key_header = ["Gene_Symbol", "Mutation_Key", "Splicing_Motif_Pos", "Motif_Type", "Is_Canonical",
              "Splicing_Key", "Splicing_Type", "Is_Inframe", "LawrenceEtAl_2014", "Cancer_Gene_Census", "Is_Cancer_Gene"]


hout = open(output_file, 'w')
print >> hout, '\t'.join(key_header) + '\t' + "Count"
with open(input_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    for i in range(len(header)):
        header2ind[header[i]] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = '\t'.join([F[header2ind[x]] for x in key_header])

        if key not in key2count: key2count[key] = 0
        key2count[key] = key2count[key] + 1


for key, value in sorted(key2count.items(), key= lambda x: x[1], reverse = True):
    if value >= 3:
        print >> hout, key + '\t' + str(value)

hout.close()                                   


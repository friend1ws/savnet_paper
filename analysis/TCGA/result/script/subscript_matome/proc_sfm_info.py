#! /usr/bin/env python

import sys, re

input_file = sys.argv[1]
output_file = sys.argv[2]

hout = open(output_file, 'w')
print >> hout, '\t'.join(["Sample name", "Cancer type", "Chr", "Start", "End", "Ref", "Alt", "Gene", "Mutation type", "Annotation"])

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        if F[8] == "U2AF1":
            pass

        F[2] = "chr" + F[2]

        annos = F[11].split(',')
        if len(annos) > 1:
            F[11] = annos[1]
    
        F[11] = re.sub(r'^\w+:', '', F[11])
        if F[10] == "": F[10] = F[7]
        if F[11] == "": F[11] = F[9]

        print >> hout, '\t'.join(F[:7]) + '\t' + F[8] + '\t' + '\t'.join(F[10:])

hout.close()


#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

hout = open(output_file, 'w')
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split('\t')
    print >> hout, '\t'.join(["Cancer type", "Sample name", "SNV count", "Indel count", "SAV count"])

    for line in hin:
        F = line.rstrip('\n').split('\t')
        print >> hout, '\t'.join(F)

hout.close()



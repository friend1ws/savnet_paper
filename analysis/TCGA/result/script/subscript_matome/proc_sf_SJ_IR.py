#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

sptype = {"exon-skip": "Exon skipping",
          "intron-retention": "Intron retention",
          "alternative-5'-splice-site": "Alternative 5'SS",
          "alternative-3'-splice-site": "Alternative 3'SS"}

hout = open(output_file, 'w')
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split('\t')
    print >> hout, '\t'.join(["Cancer type", "Gene", "Splicing position", "Splicing outcome", 
                              "Associated splicing factor mutation", "Minus Log P-value", "Q-value"])
 
    for line in hin:
        F = line.rstrip('\n').split('\t')
        print >> hout, '\t'.join([F[0], F[3], "chr" + F[2], sptype[F[4]], F[1], F[5], F[6]])

 

hout.close()


#! /usr/bin/env python

import sys

with open(sys.argv[1], 'r') as hin:
    header = hin.readline()
    for line in hin:
        F = line.rstrip('\n').split('\t')
        F[1] = str(int(F[1]) - 1)
        print '\t'.join(F)


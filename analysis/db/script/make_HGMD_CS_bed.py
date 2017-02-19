#! /usr/bin/env python

import sys, gzip

with gzip.open("/home/omega3/database/HGMD_PRO_2015.4.bed.gz", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if len(F) < 5:
            print sys.stderr, '\t'.join(F)
            continue

        if F[5].startswith("CS"): 
            print '\t'.join(F)




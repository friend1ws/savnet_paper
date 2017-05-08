#! /usr/bin/env python

import sys, gzip

input_file = sys.argv[1]
gene_symbol = sys.argv[2]


with gzip.open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[12] != gene_symbol: continue

        amino_start = int(F[6])
        amino_end = int(F[7])

        exon_starts = [int(x) for x in F[9].split(',')[:-1]]
        exon_ends = [int(x) for x in F[10].split(',')[:-1]]

        tsize = 0
        for i in range(len(exon_starts)):
            if exon_ends[i] < amino_start: continue
            if exon_starts[i] > amino_end: continue
            tsize = tsize + min(exon_ends[i], amino_end) - max(exon_starts[i], amino_start)
        
        print F[1] + '\t' + str(tsize) + '\t' + str(tsize / 3)



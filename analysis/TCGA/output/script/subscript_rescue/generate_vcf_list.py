#! /usr/bin/env python

import sys

input_file = sys.argv[1]

mut_key2gene = {}
with open(input_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        mut_key2gene[F[header2ind["Mutation_Key"]]] = F[header2ind["Gene_Symbol"]]


for mut_key in sorted(mut_key2gene):
    FF = mut_key.split(',')
    gene_symbol = mut_key2gene[mut_key]

    # deletion
    if len(FF[2]) > 1:
        print '\t'.join([FF[0], str(int(FF[1]) + 1), str(int(FF[1]) + len(FF[2]) - 1), FF[2][1:], "-"]) + '\t' + gene_symbol
    # insertion
    elif len(FF[3]) > 1:
        print '\t'.join([FF[0], FF[1], FF[1], "-", FF[3][1:]]) + '\t' + gene_symbol
    # SNV
    else:
        print '\t'.join([FF[0], FF[1], FF[1], FF[2], FF[3]]) + '\t' + gene_symbol





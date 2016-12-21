#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

key2info = {}
key2read_num = {}
with open(input_file, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["CNA"]] == "NA": continue

        key = F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Gene_Symbol"]] + '\t' + F[header2ind["Mutation_Key"]]

        is_inframe = "1" if F[header2ind["Is_Inframe"]] == "in-frame" else "0"
        is_del = "1" if int(F[header2ind["CNA"]]) < 0 else "0"

        info = is_inframe + is_del
        read_num = int(F[header2ind["Supporting_Read_Num"]])

        if key in key2info:
            if read_num > key2read_num[key]:
                key2info[key] = info
                key2read_num[key] = read_num
        else:
            key2info[key] = info
            key2read_num[key] = read_num


gene2count00 = {}
gene2count01 = {}
gene2count10 = {}
gene2count11 = {}
for key in sorted(key2info):
    sample_name, gene, mut_key = key.split('\t')

    if gene not in gene2count00:
        gene2count00[gene], gene2count01[gene], gene2count10[gene], gene2count11[gene] = 0, 0, 0, 0

    if key2info[key] == "00": gene2count00[gene] = gene2count00[gene] + 1
    if key2info[key] == "01": gene2count01[gene] = gene2count01[gene] + 1
    if key2info[key] == "10": gene2count10[gene] = gene2count10[gene] + 1
    if key2info[key] == "11": gene2count11[gene] = gene2count11[gene] + 1


hout = open(output_file, 'w')
print >> hout, "Gene_Symbol" + '\t' + "X00\tX01\tX10\tX11"
for gene in sorted(gene2count00):
    print >> hout, gene + '\t' + str(gene2count00[gene]) + '\t' + str(gene2count01[gene]) + '\t' + str(gene2count10[gene]) + '\t' + str(gene2count11[gene])

hout.close()







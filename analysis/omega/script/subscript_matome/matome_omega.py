#! /usr/bin/env python

import sys, os, glob, re

input_dir = sys.argv[1]
output_file_target = sys.argv[2]
cancer_gene_file = sys.argv[3]
BF_thres = sys.argv[4]
q_value_thres = sys.argv[5]

gene2cancer_gene_info = {}
with open(cancer_gene_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        cancer_gene_flag = "FALSE"
        if F[1] != "---" or F[2] not in ["---", "T"]:
            cancer_gene_flag = "TRUE"
        gene2cancer_gene_info[F[0]] = F[1] + '\t' + F[2] + '\t' + cancer_gene_flag



allfiles_target = glob.glob(input_dir + "/*/*.genomon_splicing_mutation.result.txt")


hout = open(output_file_target, 'w')

header_print_flag = 0
header2ind = {}
for gsfile in allfiles_target:

    cancer_type = os.path.basename(gsfile).replace(".genomon_splicing_mutation.result.txt", "")

    with open(gsfile, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        if header_print_flag == 0:
            print >> hout, "Cancer_Type" + '\t' + '\t'.join(header) + '\t' + "LawrenceEtAl_2014" + '\t' + "Cancer_Gene_Census" + '\t' + "Is_Cancer_Gene" 
            header_print_flag = 1

        for line in hin:
            F = line.rstrip('\n').split('\t')

            gene_symbol = F[0]
            sample_names = F[1].split(';')

            if float(F[header2ind["Score"]]) < float(BF_thres): continue
            if float(F[header2ind["Q_Value"]]) > float(q_value_thres): continue

            cancer_gene_info = gene2cancer_gene_info[gene_symbol] if gene_symbol in gene2cancer_gene_info else "---\t---\tFALSE" 

            for sample_name in sample_names:
                print >> hout, cancer_type + '\t' + gene_symbol + '\t' + sample_name + '\t' + '\t'.join(F[2:]) + '\t' + cancer_gene_info


hout.close()



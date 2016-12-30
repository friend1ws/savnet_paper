#! /usr/bin/env python

import sys, os, gzip, glob, re

input_dir = sys.argv[1]
output_file_target = sys.argv[2]
cancer_gene_file = sys.argv[3]
hgmd_file = sys.argv[4]
BF_thres = sys.argv[5]
q_value_thres = sys.argv[6]

gene2cancer_gene_info = {}
with open(cancer_gene_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        cancer_gene_flag = "FALSE"
        # if F[1] != "---" or F[2] not in ["---", "T"]:
        if F[4] == "CG":
            cancer_gene_flag = "TRUE"
        gene2cancer_gene_info[F[0]] = F[1] + '\t' + F[2] + '\t' + F[4] + '\t' + cancer_gene_flag

key2hgmd = {}
with gzip.open(hgmd_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        FF = F[5].split(';')
        key = F[0].replace("chr", '') + ',' + F[2] + ',' + F[3] + ',' + F[4]
        if key not in key2hgmd: key2hgmd[key] = []
        key2hgmd[key].append(FF[-1].strip('"'))


allfiles_target = glob.glob(input_dir + "/*/*.genomon_splicing_mutation.result.txt")


hout = open(output_file_target, 'w')

header_print_flag = 0
header2ind = {}
for gsfile in sorted(allfiles_target):

    # cancer_type = os.path.basename(gsfile).replace(".genomon_splicing_mutation.result.txt", "")
    cancer_type = os.path.basename(os.path.dirname(gsfile))
    if "bk" in cancer_type: continue
    if cancer_type in ["LAML", "OV"]: continue
    print gsfile

    with open(gsfile, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        if header_print_flag == 0:
            print >> hout, "Cancer_Type" + '\t' + '\t'.join(header) + '\t' + \
                           "LawrenceEtAl_2014" + '\t' + "Cancer_Gene_Census" + '\t' + "YeEtAl_2016" + '\t' + "Is_Cancer_Gene" + '\t' + "HGMD" 
            header_print_flag = 1

        for line in hin:
            F = line.rstrip('\n').split('\t')

            gene_symbol = F[0]
            sample_names = F[1].split(';')
            sp_counts = F[header2ind["Supporting_Read_Num"]].split(';')

            if float(F[header2ind["Score"]]) < float(BF_thres): continue
            if float(F[header2ind["Q_Value"]]) > float(q_value_thres): continue

            cancer_gene_info = gene2cancer_gene_info[gene_symbol] if gene_symbol in gene2cancer_gene_info else "---\t---\t---\tFALSE" 
            hgmd_info_print = ';'.join(list(set(key2hgmd[F[header2ind["Mutation_Key"]]]))) if F[header2ind["Mutation_Key"]] in key2hgmd else "---"

            for i in range(len(sample_names)):
                try:
                    print >> hout, cancer_type + '\t' + gene_symbol + '\t' + sample_names[i] + '\t' + '\t'.join(F[2:9]) + '\t' + \
                                    sp_counts[i] + '\t' + F[header2ind["Score"]] + '\t' + F[header2ind["Q_Value"]] + '\t' + \
                                    cancer_gene_info + '\t' + hgmd_info_print
                except Exception as inst:
                    print >> hout, cancer_type + '\t' + gene_symbol + '\t' + sample_names[i] + '\t' + '\t'.join(F[2:9]) + '\t' + \
                                    sp_counts[0] + '\t' + F[header2ind["Score"]] + '\t' + F[header2ind["Q_Value"]] + '\t' + \
                                    cancer_gene_info + '\t' + hgmd_info_print

hout.close()



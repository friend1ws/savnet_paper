#! /usr/bin/env python

import sys, glob, os

# expression
allfiles = glob.glob("/home/eva/genomon_out/rna_2_4_0/ogawalab/ATL/expression/*/*.sym2fkpm.txt")
key2exp = {}
sample2rna_seq = {}
for exp_file in sorted(allfiles):
    sample = os.path.basename(exp_file).replace(".sym2fkpm.txt", "")
    with open(exp_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key2exp[sample + 'T' + '\t' + F[0]] = F[1]
            sample2rna_seq[sample + 'T'] = 1


# splicing mutation
gsm_file =  "../output/ATL_EXS.genomon_splicing_mutation.result.txt"
header2ind = {}
gene_sample2score = {}
with open(gsm_file) as hin:
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene_sample2score[F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Gene_Symbol"]]] = F[header2ind["Score"]] + '\t' + F[header2ind["Q_Value"]]


allfiles = glob.glob("/home/yshira/mypaper/genomonSV_paper/analysis/ATL/result/genomon_2_3_0/Hiseq/mutation/*/*.genomon_mutation.result.filt.txt") + \
           glob.glob("/home/yshira/mypaper/genomonSV_paper/analysis/ATL/result/genomon_2_3_0/Xten_add2/mutation/*/*.genomon_mutation.result.filt.txt")


hout = open("ATL_EXS.splicing_mutation.exp_gsm.txt", 'w')
for mut_file in sorted(allfiles):
    sample = os.path.basename(os.path.dirname(mut_file))
    with open(mut_file, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')
            key = sample + '\t' +  F[6]
            if F[5] == "splicing":
                if sample in sample2rna_seq:
                    exp = round(float(key2exp[key]),4) if key in key2exp else "0.000"
                    gsm_score = gene_sample2score[key] if key in gene_sample2score else "---\t---"
                    print >> hout, sample + '\t' + '\t'.join(F[:7]) + '\t' + str(exp) + '\t' + gsm_score

hout.close()



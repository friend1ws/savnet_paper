#! /usr/bin/env python

import sys, os, glob, math

savnet_input_file = sys.argv[1]
savnet_output_file = sys.argv[2]
exp_dir = sys.argv[3]
output_file = sys.argv[4]


cancer_type = os.path.basename(savnet_input_file).replace(".mut_SJ_IR_list.txt", "")

key2func = {}
key2rank = {} # 
func2rank =  {"nonframeshift deletion": 0, \
              "nonframeshift insertion": 1, \
              "synonymous SNV": 2, \
              "nonsynonymous SNV": 3, \
              "frameshift deletion": 4, \
              "frameshift insertion": 5, \
              "stopgain": 6, \
              "splicing": 7, \
              "GSM": 8}

##########
# reading genomon splicing mutation output
key2func_gsm = {}
with open(savnet_output_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["Cancer_Type"]] != cancer_type: continue
        if F[header2ind["IR_filtered"]] == "TRUE": continue
        key = F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Gene_Symbol"]]
        if key not in key2func_gsm: 
            key2func_gsm[key] = F[header2ind["Splicing_Class"]] + '\t' + F[header2ind["Is_Inframe"]]
        else:
            key2func_gsm[key] = "complex" + '\t' + "---"
##########

##########
# reading genomon splicing mutation input for getting mutation infos
exp_barcode_list = {}
with open(savnet_input_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i


    for line in hin:

        F = line.rstrip('\n').split('\t')
        sample_name = F[header2ind["Sample_Name"]]
        mut_file = F[header2ind["Mutation_File"]]
        exp_barcode_name = os.path.basename(F[header2ind["SJ_File"]]).replace(".SJ.out.tab", "")
        exp_barcode_list[exp_barcode_name] = 1

        with open(mut_file, 'r') as hin2:

            header2ind2 = {}
            for line2 in hin2:
                if line2.startswith("#"): continue 
                if line2.startswith('\t'.join(["Chr", "Start", "End", "Ref", "Alt"])):
                    header2 = line2.rstrip('\n').split('\t')
                    for (i, cname) in enumerate(header2):
                        header2ind2[cname] = i
                    continue

                F2 = line2.rstrip('\n').split('\t')
                if "," in F2[header2ind2["Gene.refGene"]]: continue # skip mutations with multiple gene annotations"
                # if float(F2[header2ind["misRate_tumor"]]) < 0.1: continue # skip mutations with low allele freq.

                key = sample_name + '\t' + F2[header2ind2["Gene.refGene"]]

                if F2[header2ind2["Func.refGene"]] == "splicing":
                    if key not in key2func or func2rank["splicing"] > func2rank[key2func[key]]:
                        key2func[key] = "splicing"


                if F2[header2ind2["Func.refGene"]] != "exonic": continue

                if F2[header2ind2["ExonicFunc.refGene"]] in ["synonymous SNV", "nonsynonymous SNV", "stopgain", "nonframeshift deletion", "nonframeshift insertion", \
                                                            "frameshift deletion", "frameshift insertion"]:
                    if key not in key2func or func2rank[F2[header2ind2["ExonicFunc.refGene"]]] > func2rank[key2func[key]]:
                        key2func[key] = F2[header2ind2["ExonicFunc.refGene"]]
##########


all_files = glob.glob(exp_dir + "/*.sym2fpkm.txt")

gene2exp_sum = {}
gene2exp_sum2 = {}
key2fpkm = {}
sample_num = 0
for exp_file in sorted(all_files):
    sample_num = sample_num + 1
    barcode_name = os.path.basename(exp_file).replace(".sym2fpkm.txt", '')
    if barcode_name not in exp_barcode_list: continue
    sample_name = barcode_name[:15]
    # print sample_name
    with open(exp_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] not in gene2exp_sum: gene2exp_sum[F[0]] = 0.0
            if F[0] not in gene2exp_sum2: gene2exp_sum2[F[0]] = 0.0

            key = sample_name + '\t' + F[0]
            if key in key2func_gsm or key in key2func:
                key2fpkm[key] = float(F[1]) 
  
            gene2exp_sum[F[0]] = gene2exp_sum[F[0]] + float(F[1])
            gene2exp_sum2[F[0]] = gene2exp_sum2[F[0]] + float(F[1]) * float(F[1])


gene2mean = {}
gene2std = {}
for gene in sorted(gene2exp_sum):

    mean_exp = gene2exp_sum[gene] / sample_num    
    gene2mean[gene] = mean_exp
    gene2std[gene] = math.sqrt(gene2exp_sum2[gene] / sample_num - mean_exp * mean_exp)



hout = open(output_file, 'w')
for key in list(set(key2func_gsm.keys() + key2func.keys())):

    sample_name, gene = key.split('\t')
    if gene not in gene2mean: continue
    fpkm = round(key2fpkm[key], 3) if key in key2fpkm else 0.000
    if gene2std[gene] <= 0.001: continue

    fpkm_nrm = round((fpkm - gene2mean[gene]) / gene2std[gene], 3)

    if key in key2func_gsm:
        print >> hout, key + '\t' + key2func_gsm[key] + '\t' + "TRUE" + '\t' + \
                       str(fpkm) + '\t' + str(fpkm_nrm) + '\t' + str(round(gene2mean[gene], 3)) + '\t' + str(round(gene2std[gene], 3))
    else:
        print >> hout, key + '\t' + key2func[key] + '\t' + "---" + '\t' + "FALSE" + '\t' + \
                       str(fpkm) + '\t' + str(fpkm_nrm) + '\t' + str(round(gene2mean[gene], 3)) + '\t' + str(round(gene2std[gene], 3))

 
hout.close()


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
        if F[1] != "---" or F[2] != "---" or F[3] != "---" or F[4] != "---":
            cancer_gene_flag = "TRUE"
        gene2cancer_gene_info[F[0]] = F[1] + '\t' + F[2] + '\t' + F[3] + '\t' + F[4] + '\t' + cancer_gene_flag


key2hgmd = {}
with gzip.open(hgmd_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        FF = F[5].split(';')
        key = F[0].replace("chr", '') + ',' + F[2] + ',' + F[3] + ',' + F[4]
        if key not in key2hgmd: key2hgmd[key] = []
        key2hgmd[key].append(FF[-1].strip('"'))


allfiles_target = glob.glob(input_dir + "/*/*.savnet.result.txt")


"""
mut_key2mut_type = {}
mut_key2is_canonical = {}
for gsfile in sorted(allfiles_target):

    cancer_type = os.path.basename(os.path.dirname(gsfile))

    with open(gsfile, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i
 
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[header2ind["Mutation_Key"]] not in mut_key2mut_type or "creation" in mut_key2mut_type[F[header2ind["Mutation_Key"]]]:
                mut_key2mut_type[F[header2ind["Mutation_Key"]]] = F[header2ind["Mutation_Type"]]

            if F[header2ind["Mutation_Key"]] not in mut_key2is_canonical or mut_key2is_canonical[F[header2ind["Mutation_Key"]]] == "non-canonical":
                mut_key2is_canonical[F[header2ind["Mutation_Key"]]] = F[header2ind["Is_Canonical"]] 
"""

key2int_ss_by_dis = {}
key2intron = {}
for gsfile in sorted(allfiles_target):
    cancer_type = os.path.basename(os.path.dirname(gsfile))
    with open(gsfile, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')

            key = F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Mutation_Key"]]

            if F[header2ind["Splicing_Class"]] in ["alternative-5'-splice-site", "alternative-3'-splice-site", \
                                                   "intronic-alternative-5'-splice-site", "intronic-alternative-3'-splice-site"]:

                pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', F[header2ind["Splicing_Key"]])
                schr, sstart, send = pos_match.group(1), pos_match.group(2), pos_match.group(3)

                pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)\,([\+\-])', F[header2ind["Motif_Pos"]])
                mchr, mstart, mend, mdir = pos_match.group(1), pos_match.group(2), pos_match.group(3), pos_match.group(4)

                if "alternative-5'-splice-site" in F[header2ind["Splicing_Class"]] and F[header2ind["Mutation_Type"]] == "splicing donor disruption":
                    if (mdir == "+" and int(mend) < int(sstart)) or (mdir == "-" and int(send) < int(mstart)):
                        key2int_ss_by_dis[key] = 1
                        # print '\t'.join(F)
                elif "alternative-3'-splice-site" in F[header2ind["Splicing_Class"]] and F[header2ind["Mutation_Type"]] == "splicing acceptor disruption":
                    if (mdir == "-" and int(mend) < int(sstart)) or (mdir == "+" and int(send) < int(mstart)):
                        key2int_ss_by_dis[key] = 1
                        # print '\t'.join(F)


            if F[header2ind["Splicing_Class"]] == "intron-retention":
                key2intron[key] = 1 



hout = open(output_file_target, 'w')

header_print_flag = 0
for gsfile in sorted(allfiles_target):

    cancer_type = os.path.basename(os.path.dirname(gsfile))

    with open(gsfile, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        if header_print_flag == 0:
            print >> hout, "Cancer_Type" + '\t' + '\t'.join(header) + '\t' + "IR_filtered" + '\t' + \
                           "LawrenceEtAl_2014" + '\t' + "Cancer_Gene_Census" + '\t' + "VogelsteinEtAl_2013" + '\t' + "YeEtAl_2016" + '\t' + \
                           "Is_Cancer_Gene" + '\t' + "HGMD" 
            header_print_flag = 1

        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Mutation_Key"]]

            IR_filtered = "FALSE"
            if F[header2ind["Splicing_Class"]] == "intron-retention" and key in key2int_ss_by_dis: IR_filtered = "TRUE"
            if F[header2ind["Splicing_Class"]] == "opposite-side-intron-retention" and key in key2intron: IR_filtered = "TRUE"  


            gene_symbol = F[0]
            sample_names = F[1].split(';')
            sp_counts = F[header2ind["Supporting_Read_Num"]].split(';')

            # F[header2ind["Mutation_Type"]] = mut_key2mut_type[F[header2ind["Mutation_Key"]]]
            # F[header2ind["Is_Canonical"]] = mut_key2is_canonical[F[header2ind["Mutation_Key"]]]

            if float(F[header2ind["Score"]]) < float(BF_thres): continue
            if float(F[header2ind["Q_Value"]]) > float(q_value_thres): continue

            cancer_gene_info = gene2cancer_gene_info[gene_symbol] if gene_symbol in gene2cancer_gene_info else "---\t---\t---\t---\tFALSE" 
            hgmd_info_print = ';'.join(list(set(key2hgmd[F[header2ind["Mutation_Key"]]]))) if F[header2ind["Mutation_Key"]] in key2hgmd else "---"

            for i in range(len(sample_names)):
                try:
                    print >> hout, cancer_type + '\t' + gene_symbol + '\t' + sample_names[i] + '\t' + '\t'.join(F[2:9]) + '\t' + \
                                    sp_counts[i] + '\t' + F[header2ind["Score"]] + '\t' + F[header2ind["Q_Value"]] + '\t' + \
                                    IR_filtered + '\t' + cancer_gene_info + '\t' + hgmd_info_print
                except Exception as inst:
                    print >> hout, cancer_type + '\t' + gene_symbol + '\t' + sample_names[i] + '\t' + '\t'.join(F[2:9]) + '\t' + \
                                    sp_counts[0] + '\t' + F[header2ind["Score"]] + '\t' + F[header2ind["Q_Value"]] + '\t' + \
                                    IR_filtered + '\t' + cancer_gene_info + '\t' + hgmd_info_print

hout.close()




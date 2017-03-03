#! /usr/bin/env python

import sys, os, glob, re

gsm_input_dir = sys.argv[1]
gsm_output = sys.argv[2]
hotspot_result = sys.argv[3]

mut2splicing = {}
key2gsm = {}
sample_gene2gsm = {}
mut_sp2info = {}
header2ind_gsm = {}
with open(gsm_output, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    header2ind_gsm = header2ind
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["IR_filtered"]] == "TRUE": continue

        mut = F[header2ind["Mutation_Key"]]
        if mut not in mut2splicing: mut2splicing[mut] = []
        mut2splicing[mut].append(F[header2ind["Splicing_Key"]])
        
        key = F[header2ind["Sample_Name"]] + '\t' + mut + '\t' + F[header2ind["Splicing_Key"]]
        key2gsm[key] = 1

        mut_sp2info[mut + '\t' + F[header2ind["Splicing_Key"]]] = '\t'.join(F)

        sample_gene = F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Gene_Symbol"]]
        sample_gene2gsm[sample_gene] = 1


for mut in mut2splicing:
    splicing = list(set(mut2splicing[mut]))
    mut2splicing[mut] = splicing

sample2SJ_file = {}
sample2IR_file = {}
gsm_input_files = glob.glob(gsm_input_dir + "/*.mut_SJ_IR_list.txt")
sample_mut2exists = {}
for gsm_input in sorted(gsm_input_files):
    
    cancer_type = os.path.basename(gsm_input).replace(".mut_SJ_IR_list.txt", "")

    with open(gsm_input, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            sample = F[header2ind["Sample_Name"]]
            sample2SJ_file[sample] = F[header2ind["SJ_File"]]
            sample2IR_file[sample] = F[header2ind["IR_File"]]

            mut_file = F[header2ind["Mutation_File"]]
            print >> sys.stderr, "Reading mutation file: " + cancer_type + '\t' + sample
            with open(mut_file, 'r') as hin2:
                header2ind2 = {}
                header2 = hin2.readline().rstrip('\n').split('\t')
                for (i, cname) in enumerate(header2):
                   header2ind2[cname] = i

                for line2 in hin2:
                    F2 = line2.rstrip('\n').split('\t')
                    mut = ','.join([F2[header2ind2[x]] for x in ["Chr", "End", "Ref", "Alt"]])
                    sample_mut2exists[sample + '\t' + mut] = 1


with open(hotspot_result, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        sample = F[header2ind["Sample_Name"]]
        mut = ','.join([F[header2ind[x]] for x in ["Chr", "End", "Ref", "Alt"]])
        cancer_type = F[header2ind["Cancer_Type"]]
        gene = F[header2ind["Gene_Symbol"]]
        # pvalue = F[header2ind["GenomonPval(EBCall)"]]

        if sample + '\t' + mut in sample_mut2exists: 
            continue
 
        if sample not in sample2SJ_file: continue
        if sample not in sample2IR_file: continue
        if mut not in mut2splicing: continue
        for splicing in mut2splicing[mut]:
            # print sample + '\t' + mut + '\t' + splicing

            key = sample + '\t' + mut + '\t' + splicing
            is_gsm_key = "TRUE" if key in key2gsm else "FALSE"
            is_gsm_gene = "TRUE" if sample + '\t' + gene in sample_gene2gsm else "FALSE"

            if is_gsm_key == "TRUE": continue # sample + mut + splicing level check
            if is_gsm_gene == "TRUE": continue # sample + gene level check

            sp_count = 0
            pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', splicing)
            schr, sstart, send = pos_match.group(1), pos_match.group(2), pos_match.group(3)

            # SJ file
            if sstart != send:
                with open(sample2SJ_file[sample], 'r') as hin2:
                    for line2 in hin2:
                        F2 = line2.rstrip('\n').split('\t')
                        if schr == F2[0] and sstart == F2[1] and send == F2[2]:
                            sp_count = int(F2[6])
                     
            # IR_file   
            if sstart == send:
                with open(sample2IR_file[sample], 'r') as hin2:
                    header2ind2 = {}
                    header2 = hin2.readline().rstrip('\n').split('\t')
                    for (i, cname) in enumerate(header2):
                        header2ind2[cname] = i
        
                    for line2 in hin2:
                        F2 = line2.rstrip('\n').split('\t')
                        if schr == F2[header2ind2["Chr"]] and sstart == F2[header2ind2["Boundary_Pos"]]:
                            sp_count = int(F2[header2ind2["Intron_Retention_Read_Count"]])


            if sp_count <= 2: continue

            infos = mut_sp2info[mut + '\t' + splicing].split('\t')

            if sp_count > 0:
                print cancer_type + '\t' + gene + '\t' + sample + '\t' + \
                       '\t'.join([infos[header2ind_gsm[x]] for x in ["Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical", 
                                                                     "Splicing_Key", "Splicing_Class", "Is_Inframe"]]) + '\t' + \
                       str(sp_count) + '\t' + "Rescued" + '\t' + "Rescued" + '\t' + \
                       '\t'.join([infos[header2ind_gsm[x]] for x in ["IR_filtered", "LawrenceEtAl_2014", "Cancer_Gene_Census", "VogelsteinEtAl_2013", "YeEtAl_2016", 
                                                                     "Is_Cancer_Gene", "HGMD"]])


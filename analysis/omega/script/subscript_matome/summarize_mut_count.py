#! /usr/bin/env python

import sys, os, glob

gsm_file = sys.argv[1]
input_list_dir = sys.argv[2]
output_file = sys.argv[3]

gsm_key = {}
with open(gsm_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        gsm_key[F[header2ind["Cancer_Type"]] + '\t' + F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Mutation_Key"]]] = 1


sample2gsm_count = {}
for gsm_key in sorted(gsm_key):
    cancer_type, sample_name, mut_key = gsm_key.split('\t')
    if cancer_type + '\t' + sample_name not in sample2gsm_count: sample2gsm_count[cancer_type + '\t' + sample_name] = 0
    sample2gsm_count[cancer_type + '\t' + sample_name] = sample2gsm_count[cancer_type + '\t' + sample_name] + 1



all_list_files = glob.glob(input_list_dir + "/*.mut_SJ_IR_list.txt")
hout = open(output_file, 'w')

print >> hout, "Cancer_Type" + '\t' + "Sample_Name" + '\t' + "SNV_Count" + '\t' + "Indel_Count" + '\t' + "SASM_Count"
for list_file in sorted(all_list_files):

    cancer_type = os.path.basename(list_file).replace(".mut_SJ_IR_list.txt", "")
    if cancer_type in ["LAML", "OV"]: continue
    with open(list_file, 'r') as hin1:

        mut_file_ind = -1
        header = hin1.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            if cname == "Mutation_File": mut_file_ind = i


        for line1 in hin1:
            F1 = line1.rstrip('\n').split('\t')

            snv_count = 0
            indel_count = 0
            # mut_count = 0
            with open(F1[mut_file_ind], 'r') as hin2:
                for line2 in hin2:
                    F2 = line2.rstrip('\n').split('\t')
                    if F2[0].startswith('#'): continue
                    if F2[0] == "Chr": continue
                    if F2[3] != "-" and F2[4] != "-":
                        snv_count = snv_count + 1
                    else:
                        indel_count = indel_count + 1 
                    # mut_count = mut_count + 1

            gsm_count = sample2gsm_count[cancer_type + '\t' + F1[0]] if cancer_type + '\t' + F1[0] in sample2gsm_count else 0

            print >> hout, cancer_type + '\t' + F1[0] + '\t'+ str(snv_count) + '\t' + str(indel_count) + '\t' + str(gsm_count)

                
hout.close()




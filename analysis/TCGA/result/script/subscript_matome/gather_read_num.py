#! /usr/bin/env python

import sys, re, glob

gsm_file = sys.argv[1]
gsm_out_dir = sys.argv[2]
output_file = sys.argv[3]
motif_count_thres = 5 

sample_motif = {}
motif2splicing_info = {}
with open(gsm_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        sample_motif[F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Motif_Pos"]]] = 1

        if F[header2ind["Motif_Pos"]] not in motif2splicing_info: motif2splicing_info[F[header2ind["Motif_Pos"]]] = []
    
        splicing_info = '\t'.join(F[header2ind[x]] for x in ["Splicing_Key", "Splicing_Class", "Is_Inframe"])
 
        if splicing_info not in motif2splicing_info[F[header2ind["Motif_Pos"]]]:
            motif2splicing_info[F[header2ind["Motif_Pos"]]].append(splicing_info)


# count motif count
motif2count = {}
for sm in sample_motif:
    sample, motif = sm.split('\t')
    if motif not in motif2count: motif2count[motif] = 0
    motif2count[motif] = motif2count[motif] + 1


# check splicing files
sample2SJ_file = {}
sample2IR_file = {}
sample2weight = {}
gsm_input_files = glob.glob(gsm_out_dir + "/*.mut_SJ_IR_list.txt")

for gsm_input in sorted(gsm_input_files):
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
            sample2weight[sample] = F[header2ind["Weight"]]


processed_keys = {}
hout = open(output_file, 'w')
print >> hout, '\t'.join(["Cancer_Type", "Gene_Symbol", "Sample_Name", "Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical",
                          "Splicing_Key", "Splicing_Class", "Is_Inframe", "Supporting_Read_Num", "Weight"])

with open(gsm_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if int(motif2count[F[header2ind["Motif_Pos"]]]) <= motif_count_thres: continue

        temp_key = '\t'.join(F[header2ind[x]] for x in ["Cancer_Type", "Gene_Symbol", "Sample_Name", 
                            "Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical"])

        if temp_key in processed_keys: continue
        processed_keys[temp_key] = 1

        sample = F[header2ind["Sample_Name"]]
        motif_pos = F[header2ind["Motif_Pos"]]
        weight = sample2weight[sample]

        for sp_info in motif2splicing_info[motif_pos]:
            sp_key, sp_class, is_inframe = sp_info.split('\t')
            pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', sp_key)
            schr, sstart, send = pos_match.group(1), pos_match.group(2), pos_match.group(3)

            sp_count = 0
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

            print >> hout, temp_key + '\t' + sp_key + '\t' + sp_class + '\t' + is_inframe + '\t' + str(sp_count) + '\t' + weight



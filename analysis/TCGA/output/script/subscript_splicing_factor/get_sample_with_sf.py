#! /usr/bin/env python

import sys, glob, os

# mut_dir = sys.argv[1]
# all_file = glob.glob(mut_dir + "/*/*.mutation.filt.txt")

input_list_dir = sys.argv[1]
sfgene_num_thres = int(sys.argv[2])
sfgene_ratio_thres = float(sys.argv[3])
sfgene_list_output = sys.argv[4]
cancer_type_sfgene_list = sys.argv[5]


all_files = glob.glob(input_list_dir + "/*.mut_SJ_IR_list.txt")

hout1 = open(sfgene_list_output, 'w')
hout2 = open(cancer_type_sfgene_list, 'w')
for input_list_file in sorted(all_files):

    cancer_type = os.path.basename(input_list_file).replace(".mut_SJ_IR_list.txt", '')
    sample_count = 0

    sfgene2count = {"U2AF1": 0, "SRSF2": 0, "SF3B1": 0, "ZRSR2": 0}
    sfgene2line = {"U2AF1": [], "SRSF2": [], "SF3B1": [], "ZRSR2": []}
    with open(input_list_file) as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            sample = F[header2ind["Sample_Name"]] 
            mut_file = F[header2ind["Mutation_File"]]
            sample_count = sample_count + 1

            with open(mut_file, 'r') as hin2:
                for line2 in hin2:
                    F2 = line2.rstrip('\n').split('\t')
                    if F2[0].startswith("#"): continue
                    if F2[6] == "U2AF1" and (("S34" in F2[9]) or ("Q157" in F2[9])):
                        sfgene2count["U2AF1"] = sfgene2count["U2AF1"] + 1
                        sfgene2line["U2AF1"].append(sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10]))
                        # print sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10])

                    if F2[6] == "ZRSR2" and F2[5] in ["exonic", "splicing"] and F2[8] != "synonymous SNV":
                        sfgene2count["ZRSR2"] = sfgene2count["ZRSR2"] + 1
                        sfgene2line["ZRSR2"].append(sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10]))
                        # print sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10])

                    if F2[6] == "SF3B1" and (("K700" in F2[9]) or ("K666" in F2[9]) or ("H662" in F2[9]) or \
                                            ("R625" in F2[9]) or ("E622" in F2[9]) or ("G740" in F2[9]) or \
                                            ("G742" in F2[9]) or ("N626" in F2[9]) or ("E902" in F2[9]) or ("R957" in F2[9])):
                        sfgene2count["SF3B1"] = sfgene2count["SF3B1"] + 1
                        sfgene2line["SF3B1"].append(sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10]))
                        # print sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10])    

                    if F2[6] == "SRSF2" and "P95" in F2[9]:
                        sfgene2count["SRSF2"] = sfgene2count["SRSF2"] + 1
                        sfgene2line["SRSF2"].append(sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10]))
                        # print sample + '\t' + cancer_type + '\t' + '\t'.join(F2[:10])


    for sfgene in ["U2AF1", "SRSF2", "SF3B1", "ZRSR2"]:
        if sfgene2count[sfgene] > sfgene_num_thres and float(sfgene2count[sfgene]) / sample_count > sfgene_ratio_thres:
            print >> hout1, '\n'.join(sfgene2line[sfgene])
            print >> hout2, cancer_type + '\t' + sfgene


hout1.close()
hout2.close()


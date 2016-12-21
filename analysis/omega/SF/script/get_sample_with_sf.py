#! /usr/bin/env python

import sys, glob, os

# all_file = glob.glob("/home/kchiba/work_directory/work_hotspot/black_list_output_min10/*/*.genomon.genomon_mutation.result.filt.blacklist_filtered.txt")
all_file = glob.glob("../../mutation/output/*/*.mutation.filt.txt")

for mut_file in sorted(all_file):
    sample = os.path.basename(mut_file).replace(".mutation.filt.txt", "")
    cancer_type = os.path.basename(os.path.dirname(mut_file))
    with open(mut_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0].startswith("#"): continue
            if F[6] == "U2AF1" and (("S34" in F[9]) or ("Q157" in F[9])):
                print sample + '\t' + cancer_type + '\t' + '\t'.join(F[:10])

            if F[6] == "ZRSR2" and F[5] in ["exonic", "splicing"] and F[8] != "synonymous SNV":
                print sample + '\t' + cancer_type + '\t' + '\t'.join(F[:10])

            if F[6] == "SF3B1" and (("K700" in F[9]) or ("K666" in F[9]) or ("H662" in F[9]) or \
                                    ("R625" in F[9]) or ("E622" in F[9]) or ("G740" in F[9]) or \
                                    ("G742" in F[9]) or ("N626" in F[9]) or ("E902" in F[9]) or ("R957" in F[9])):
                print sample + '\t' + cancer_type + '\t' + '\t'.join(F[:10])    

            if F[6] == "SRSF2" and "P95" in F[9]:
                print sample + '\t' + cancer_type + '\t' + '\t'.join(F[:10])


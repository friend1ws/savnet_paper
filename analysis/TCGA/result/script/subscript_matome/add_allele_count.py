#! /usr/bin/env python

import sys, os, glob

input_file = sys.argv[1]
allele_count_dir = sys.argv[2]
output_file = sys.argv[3]

print_header = ""
header2ind = {}
key2info = {}
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split('\t')
    print_header = '\t'.join(header)
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["Splicing_Class"]] != "intron-retention": continue
        FF = F[header2ind["Mutation_Key"]].split(',')
        key = F[header2ind["Sample_Name"]] + '\t' + FF[0] + '\t' + FF[1]
        key2info[key] = '\t'.join(F)


hout = open(output_file, 'w')
print >> hout, print_header + '\t' + "Splice_Junction_Negative" + '\t' + "Splice_Junction_Positive" + '\t' + "Intron_Retention_Negative" + '\t' + "Intron_Retention_Positive"

all_files = glob.glob(allele_count_dir + "/*/*.allele_count.txt")
for allele_count_file in sorted(all_files):

    sample_name = os.path.basename(allele_count_file).replace(".allele_count.txt", '')
    header2ind2 = {}
    with open(allele_count_file) as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind2[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')

            start_mut = F[header2ind2["Start_Mut"]] if F[header2ind2["Ref_Mut"]] != "-" else str(int(F[header2ind2["Start_Mut"]]) - 1)
            key = sample_name + '\t' + F[header2ind2["Chr_Mut"]] + '\t' + start_mut
    
            if key in key2info:
                print >> hout, key2info[key] + '\t' + '\t'.join([F[header2ind2[x]] for x in ["Splice_Junction_Negative", "Splice_Junction_Positive", "Intron_Retention_Negative", "Intron_Retention_Positive"]])


hout.close()



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
        # if F[header2ind["Splicing_Class"]] != "intron-retention": continue
        FF = F[header2ind["Mutation_Key"]].split(',')
        key = F[header2ind["Sample_Name"]] + '\t' + FF[0] + '\t' + FF[1]
        key2info[key] = key2info[key] + ';' + F[header2ind["Splicing_Class"]] if key in key2info else F[header2ind["Splicing_Class"]]
        # key2info[key] = '\t'.join(F)


hout = open(output_file, 'w')
print_header = ["Gene_Symbol", "Chr_Mut", "Start_Mut", "End_Mut", "Ref_Mut", "Alt_Mut", "Chr_Motif", "Start_Motif", "End_Motif", \
          "Type_Motif", "Strand_Motif"]

print >> hout, "Cancer_Type" + '\t' + "Sample_Name" + '\t' + \
               '\t'.join(print_header + ["GenomonSplicingMutation", "Rel_Start_Motif", "Rel_End_Motif"])

# ERBB4   2       212484000       212484000       C       A       2       212484000       212484009       acceptor        -
# print >> hout, print_header + '\t' + "Splice_Junction_Negative" + '\t' + "Splice_Junction_Positive" + '\t' + "Intron_Retention_Negative" + '\t' + "Intron_Retention_Positive"

all_files = glob.glob(allele_count_dir + "*/*.allele_count.txt")
for allele_count_file in sorted(all_files):

    cancer_type = os.path.basename(os.path.dirname(allele_count_file))
    sample_name = os.path.basename(allele_count_file).replace(".allele_count.txt", '')
    header2ind2 = {}
    with open(allele_count_file) as hin:

        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind2[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[2] == "100912669":
                pass

            if F[header2ind2["Strand_Motif"]] == "+":
                rel_start = int(F[header2ind2["Start_Mut"]]) - int(F[header2ind2["Start_Motif"]]) + 1
                rel_end = int(F[header2ind2["End_Mut"]]) - int(F[header2ind2["Start_Motif"]]) + 1
            else:
                rel_start = int(F[header2ind2["End_Motif"]]) - int(F[header2ind2["End_Mut"]]) + 1
                rel_end = int(F[header2ind2["End_Motif"]]) - int(F[header2ind2["Start_Mut"]]) + 1
            rel_start = max(0, rel_start)
            rel_end = min(9, rel_end) if F[header2ind2["Type_Motif"]] == "donor" else min(10, rel_end)

            start_mut = F[header2ind2["Start_Mut"]] if F[header2ind2["Alt_Mut"]] != "-" else str(int(F[header2ind2["Start_Mut"]]) - 1)
            key = sample_name + '\t' + F[header2ind2["Chr_Mut"]] + '\t' + start_mut
   
            GSM_info = key2info[key] if key in key2info else "---" 
            print >> hout, cancer_type + '\t' + sample_name + '\t' + \
                    '\t'.join(F[header2ind2[x]] for x in print_header) + '\t' + GSM_info + '\t' + str(rel_start) + '\t' + str(rel_end)

            # if key in key2info:
            #     print >> hout, key2info[key] + '\t' + '\t'.join([F[header2ind2[x]] for x in ["Splice_Junction_Negative", "Splice_Junction_Positive", "Intron_Retention_Negative", "Intron_Retention_Positive"]])


hout.close()



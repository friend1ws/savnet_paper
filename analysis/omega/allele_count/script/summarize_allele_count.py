#! /usr/bin/env python

import sys, os, glob
import my_seq

input_file = sys.argv[1]
allele_count_dir = sys.argv[2]
reference = sys.argv[3]
output_file = sys.argv[4]

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


# get sample-gene info for getting fkpm values
sample_gene2fkpm = {}
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
            sample_gene2fkpm[sample_name + '\t' + F[header2ind2["Gene_Symbol"]]] = 0.000


# get expression info
all_files = glob.glob("/home/yshira/CD274/CD274_proj/result/TCGA_160316/*/expression/*/*.sym2fkpm.txt")
sample_isrna = {}
for exp_file in sorted(all_files):
    cancer_type = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(exp_file))))
    sample_name = os.path.basename(exp_file).replace(".sym2fkpm.txt", '')[:15]
    sample_isrna[sample_name] = 1

    with open(exp_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            sample_gene = sample_name + '\t' + F[0]
            if sample_gene in sample_gene2fkpm:
                # print cancer_type + '\t' + sample_gene
                sample_gene2fkpm[sample_gene] = round(float(F[1]), 3)


hout = open(output_file, 'w')
print_header = ["Gene_Symbol", "Chr_Mut", "Start_Mut", "End_Mut", "Ref_Mut", "Alt_Mut", "Chr_Motif", "Start_Motif", "End_Motif", \
          "Type_Motif", "Strand_Motif"]

print >> hout, "Cancer_Type" + '\t' + "Sample_Name" + '\t' + \
               '\t'.join(print_header + ["GenomonSplicingMutation", "FKPM", "Rel_Start_Motif", "Rel_End_Motif", "Motif_Seq"])

all_files = glob.glob(allele_count_dir + "*/*.allele_count.txt")
for allele_count_file in sorted(all_files):

    cancer_type = os.path.basename(os.path.dirname(allele_count_file))
    sample_name = os.path.basename(allele_count_file).replace(".allele_count.txt", '')
    if sample_name not in sample_isrna:
        print sample_name
        continue

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
                # if F[header2ind2["Type_Motif"]] == "donor":
                #     motif_start = int(F[header2ind2["Start_Motif"]]) - 1
                #     motif_end = int(F[header2ind2["End_Motif"]]) + 0
                # else:
                #     motif_start = int(F[header2ind2["Start_Motif"]]) - 12
                #     motif_end = int(F[header2ind2["End_Motif"]]) + 2
                motif_start = int(F[header2ind2["Start_Motif"]])
                motif_end = int(F[header2ind2["End_Motif"]])
                motif_pos = F[header2ind2["Chr_Motif"]] + ':' + str(motif_start) + '-' + str(motif_end)
                motif_seq = my_seq.get_seq(reference, motif_pos)
            else:
                rel_start = int(F[header2ind2["End_Motif"]]) - int(F[header2ind2["End_Mut"]]) + 1
                rel_end = int(F[header2ind2["End_Motif"]]) - int(F[header2ind2["Start_Mut"]]) + 1
                # if F[header2ind2["Type_Motif"]] == "donor":
                #     motif_start = int(F[header2ind2["Start_Motif"]]) - 0
                #     motif_end = int(F[header2ind2["End_Motif"]]) + 1
                # else:
                #     motif_start = int(F[header2ind2["Start_Motif"]]) - 2
                #     motif_end = int(F[header2ind2["End_Motif"]]) + 12
                motif_start = int(F[header2ind2["Start_Motif"]])
                motif_end = int(F[header2ind2["End_Motif"]])
                motif_pos = F[header2ind2["Chr_Motif"]] + ':' + str(motif_start) + '-' + str(motif_end)
                motif_seq = my_seq.reverse_complement(my_seq.get_seq(reference, motif_pos))

            rel_start = max(1, rel_start)
            rel_end = min(9, rel_end) if F[header2ind2["Type_Motif"]] == "donor" else min(7, rel_end)

            start_mut = F[header2ind2["Start_Mut"]] if F[header2ind2["Alt_Mut"]] != "-" else str(int(F[header2ind2["Start_Mut"]]) - 1)
            key = sample_name + '\t' + F[header2ind2["Chr_Mut"]] + '\t' + start_mut

            #        sample_gene2fkpm[sample_name + '\t' + F[header2ind2["Gene_Symbol"]]]

            fkpm = sample_gene2fkpm[sample_name + '\t' + F[header2ind2["Gene_Symbol"]]]
            # print sample_name + '\t' + F[header2ind2["Gene_Symbol"]]
            # print fkpm

            GSM_info = key2info[key] if key in key2info else "---" 
            print >> hout, cancer_type + '\t' + sample_name + '\t' + \
                    '\t'.join(F[header2ind2[x]] for x in print_header) + '\t' + GSM_info + '\t' + str(fkpm) + '\t' + \
                    str(rel_start) + '\t' + str(rel_end) + '\t' + motif_seq

            # if key in key2info:
            #     print >> hout, key2info[key] + '\t' + '\t'.join([F[header2ind2[x]] for x in ["Splice_Junction_Negative", "Splice_Junction_Positive", "Intron_Retention_Negative", "Intron_Retention_Positive"]])


hout.close()



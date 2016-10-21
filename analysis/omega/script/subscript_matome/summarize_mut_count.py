#! /usr/bin/env python

import sys, os, glob

allfiles = glob.glob("../output/*/mut_SJ_IR_list.txt")

output_file = sys.argv[1]
hout = open(output_file, 'w')

print >> hout, "Cancer_Type" + '\t' + "Sample_Name" + '\t' + "Mut_Count"
for list_file in sorted(allfiles):

    cancer_type = os.path.basename(os.path.dirname(list_file))

    with open(list_file, 'r') as hin1:

        mut_file_ind = -1
        header = hin1.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            if cname == "Mutation_File": mut_file_ind = i


        for line1 in hin1:
            F1 = line1.rstrip('\n').split('\t')
            mut_count = 0
            with open(F1[mut_file_ind], 'r') as hin2:
                for line2 in hin2:
                    F2 = line2.rstrip('\n').split('\t')
                    if F2[0].startswith('#'): continue
                    if F2[0] == "Chr": continue
                    mut_count = mut_count + 1

            print >> hout, cancer_type + '\t' + F1[0] + '\t'+ str(mut_count)

                
hout.close()




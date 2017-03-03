#! /usr/bin/env python

import sys, glob, os


# read genomon splicing mutation file
gsm_file = sys.argv[1]
input_list_dir = sys.argv[2]

header2ind0 = {}
# key2exist = {}
key2gsm = {}

with open(gsm_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind0[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind0["IR_filtered"]] == "TRUE": continue

        FF = F[header2ind0["Mutation_Key"]].split(',')
        mut_chr, mut_start, mut_end, mut_ref, mut_alt = FF[0], FF[1], FF[1], FF[2], FF[3]
        if len(mut_ref) > 1:
            mut_ref = mut_ref[1:]
            mut_alt = "-"
            mut_start = str(int(mut_start) + 1)
            mut_end = str(int(mut_start) + len(mut_ref) - 1)
        elif len(mut_alt) > 1:
            mut_ref = "-"
            mut_alt = mut_alt[1:]
        
        key = F[header2ind0["Sample_Name"]] + '\t' + mut_chr + '\t' + mut_start + '\t' + mut_end + '\t' + mut_ref + '\t' + mut_alt
        # if mut_alt == "-":
        #     print >> sys.stderr, key
        # key2exist[key] = 1

        gsm = F[header2ind0["Splicing_Class"]]
        if gsm == "intronic-alternative-5'-splice-site": gsm = "alternative-5'-splice-site"
        if gsm == "intronic-alternative-3'-splice-site": gsm = "alternative-3'-splice-site"
        if gsm == "opposite-side-intron-retention": gsm = "intron-retention"

        if key in key2gsm:
            if gsm not in key2gsm[key]:
                key2gsm[key].append(gsm)
        else:
            key2gsm[key] = [gsm]



all_list_files = glob.glob(input_list_dir + "/*.mut_SJ_IR_list.txt")

header_flag = 0
for list_file in sorted(all_list_files):

    # cancer_type = os.path.basename(os.path.dirname(list_file))
    cancer_type = os.path.basename(list_file).replace(".mut_SJ_IR_list.txt", "")
    if cancer_type == "OV" or cancer_type == "LAML": continue
    print >> sys.stderr, cancer_type
    header2ind1 = {}    
    with open(list_file, 'r') as hin1:

        header1 = hin1.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header1):
            header2ind1[cname] = i

        for line1 in hin1:
            F1 = line1.rstrip('\n').split('\t')
            sample_name = F1[header2ind1["Sample_Name"]]
            mut_file = F1[header2ind1["Mutation_File"]]

            # if sample_name in ["TCGA-BC-A112-01", "TCGA-DD-A1EG-01", "TCGA-DD-A39Y-01", "TCGA-G3-A3CJ-01"]: continue 
            header2ind2 = {}
            with open(mut_file, 'r') as hin2:

                header2 = ["#"]
                while header2[0].startswith("#"):
                    header2 = hin2.readline().rstrip('\n').split('\t')

                for (i, cname) in enumerate(header2):
                    header2ind2[cname] = i

                if header_flag == 0:
                    print "Sample_Name" + '\t' + "Cancer_Type" + '\t' + '\t'.join(header2[:7]) + '\t' + header2[8] + '\t' + "GSM"
                    header_flag = 1

                for line2 in hin2:
            
                    F2 = line2.rstrip('\n').split('\t')

                    key = sample_name + '\t' + '\t'.join(F2[header2ind2[x]] for x in ["Chr", "Start", "End", "Ref", "Alt"])
                    gsm_info = "no-change"
                    if key in key2gsm:
                        if len(key2gsm[key]) > 1:
                            gsm_info = "complex"
                        else:
                            gsm_info = key2gsm[key][0]

                    #     print >> sys.stderr, key

                    if gsm_info == "FALSE" and F2[header2ind2["Func.refGene"]] not in ["exonic", "exonic;splicing", "splicing"]: continue
                    # if gsm_info == "FALSE" and F2[header2ind2["ExonicFunc.refGene"]] == "synonymous SNV": continue

                    if F2[8].strip() == "": F2[8] = "---"

                    print sample_name + '\t' + cancer_type + '\t' + '\t'.join(F2[:7]) + '\t' + F2[8] + '\t' + gsm_info



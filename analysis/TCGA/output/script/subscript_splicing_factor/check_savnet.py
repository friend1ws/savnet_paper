#! /usr/bin/env python

import sys, os, glob

savnet_result = sys.argv[1]

# all_black_list = glob.glob("*_filt.sig.txt")

sample2sf_mut = {}
with open("../splicing_factor/sample_list_with_sfm.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] not in sample2sf_mut: sample2sf_mut[F[0]] = []
        sample2sf_mut[F[0]].append(F[8])

        
key2ctype = {}
# for bfile in sorted(all_black_list):
#     ctype = os.path.basename(bfile).replace("_SJ_IR_filt.sig.txt", "")
#     with open(bfile, 'r') as hin:
#         header = hin.readline()
#         for line in hin:
#             F = line.rstrip('\n').split('\t')
#             key2ctype[F[1].strip('"')] = ctype

SJ_black_list = glob.glob("../splicing_factor/*_SJ.sig.txt")
for bfile in sorted(SJ_black_list):
    cancer_type, sf_gene = os.path.basename(bfile).replace("_SJ.sig.txt", "").split("_")
    ctype = sf_gene + '(' + cancer_type + ')'
    with open(bfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0].strip('"') + ':' + F[1].strip('"') + '-' + F[2].strip('"')
            if key not in key2ctype: key2ctype[key] = []
            key2ctype[key].append(ctype)

IR_black_list = glob.glob("../splicing_factor/*_IR.sig.txt")
for bfile in sorted(IR_black_list):
    cancer_type, sf_gene = os.path.basename(bfile).replace("_IR.sig.txt", "").split("_")
    ctype = sf_gene + '(' + cancer_type + ')'
    with open(bfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0].strip('"') + ':' + F[1].strip('"') + '-' + F[1].strip('"')
            if key not in key2ctype: key2ctype[key] = []
            key2ctype[key].append(ctype)
            
# for key in sorted(key2ctype):
#     print key + '\t' + str(key2ctype[key])

with open(savnet_result, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    print '\t'.join(header[:9]) + '\t' + "Controled_By" + '\t' + "SF_mut_status"

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["Splicing_Key"]] in key2ctype:

            # print key2ctype[F[header2ind["Splicing_Key"]]]

            sf_mut_status = ';'.join(sample2sf_mut[F[header2ind["Sample_Name"]]]) if F[header2ind["Sample_Name"]] in sample2sf_mut else "---"
            print '\t'.join(F[:9]) + '\t' + ','.join(key2ctype[F[header2ind["Splicing_Key"]]]) + '\t' + sf_mut_status




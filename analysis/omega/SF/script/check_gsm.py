#! /usr/bin/env python

import sys, os, glob

# all_black_list = glob.glob("*_filt.sig.txt")

sample2sf_mut = {}
with open("../output/sample_list_with_sfm.txt", 'r') as hin:
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

SJ_black_list = glob.glob("../output/*_SJ.sig.txt")
for bfile in sorted(SJ_black_list):
    ctype = os.path.basename(bfile).replace(".sig.txt", "")
    with open(bfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0].strip('"') + ':' + F[1].strip('"') + '-' + F[2].strip('"')
            if key not in key2ctype: key2ctype[key] = []
            key2ctype[key].append(ctype)

IR_black_list = glob.glob("../output/*_IR.sig.txt")
for bfile in sorted(SJ_black_list):
    ctype = os.path.basename(bfile).replace(".sig.txt", "")
    with open(bfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0].strip('"') + ':' + F[1].strip('"') + '-' + F[1].strip('"')
            if key not in key2ctype: key2ctype[key] = []
            key2ctype[key].append(ctype)
            
# for key in sorted(key2ctype):
#     print key + '\t' + str(key2ctype[key])

with open("../../matome/omega.genomon_splicing_mutation.result.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[7] in key2ctype:
            sf_mut_status = ';'.join(sample2sf_mut[F[2]]) if F[2] in sample2sf_mut else "---"
            print '\t'.join(F) + '\t' + ';'.join(key2ctype[F[7]]) + '\t' + sf_mut_status



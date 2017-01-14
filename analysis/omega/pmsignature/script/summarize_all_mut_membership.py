#! /usr/bin/env python

import sys, os, glob

all_mem_files = glob.glob("../output/*/mapping/mut_membership.txt")

print '\t'.join(["Cancer_Type", "Sig_Num", "Membership_Sum", "COSM_ID", "Corr"])
for mem_file in sorted(all_mem_files):

    cosm_file = os.path.dirname(mem_file) + "/cosm_sig.txt"
    ctype = os.path.basename(os.path.dirname(os.path.dirname(mem_file)))

    sig2cosm = {}
    sig2value = {}
    with open(cosm_file, 'r') as hin:
        header = hin.readline()
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] not in sig2cosm or float(F[2]) > sig2value[F[0]]:
                sig2cosm[F[0]] = F[1]
                sig2value[F[0]] = float(F[2])


    sig2mem_sum = {}
    with open(mem_file, 'r') as hin:

        header2ind = {}
        sig_num = 0
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i
            if cname.startswith("signature_"):
                sig_num = sig_num + 1
                sig2mem_sum[str(sig_num)] = 0

        for line in hin:
            F = line.rstrip('\n').split('\t')
            for i in range(1, sig_num + 1):
                sig2mem_sum[str(i)] = sig2mem_sum[str(i)] + float(F[header2ind["signature_" + str(i)]])      
                     

    for i in sorted(sig2mem_sum):
        cosm = sig2cosm[i] if i in sig2cosm else "NA"
        cor_value = sig2value[i] if i in sig2value else "NA"
        print ctype + '\t' + str(i) + '\t' + str(sig2mem_sum[i]) + '\t' + cosm + '\t' + str(cor_value)




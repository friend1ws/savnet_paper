#! /usr/bin/env python

import sys, os, glob
import annot_utils.boundary
import subprocess


output_file = sys.argv[1]

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


omega_result_file = "../../matome/omega.motif_summary.txt"
key2gsm_info = {}
with open(omega_result_file, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["Mutation_Type"]] not in ["splicing donor creation", "splicing acceptor creation"]: continue
        motif = "donor" if F[header2ind["Mutation_Type"]] == "splicing donor creation" else "acceptor"
        key = F[header2ind["Sample_Name"]] + '\t' + "chr" + F[header2ind["Mutation_Key"]]
        key2gsm_info[key] = motif + '\t' + F[header2ind["Rel_Pos"]]
         

all_mem_files = glob.glob("../output/*/mapping/mut_membership.txt")

hout_master = open(output_file, 'w')
print >> hout_master, "Cancer_Type" + '\t' + "Splice_Site" + '\t' + "Splice_Pos" + '\t' + "Is_GSM" + '\t' + "Sig_Num" + '\t' + "Membership_Sum" + '\t' + "COSM_ID" + '\t' + "Corr"

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


    splice_pos_gsm_sig2mem_sum = {}
    with open(mem_file, 'r') as hin:

        header2ind = {}
        sig_num = 0
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i
            if cname.startswith("signature_"):
                sig_num = sig_num + 1

        for line in hin:
            F = line.rstrip('\n').split('\t')
            gsm_key1 = F[0] + '\t' + ','.join(F[1:5])
            F[3] = reverse_complement(F[3])
            F[4] = reverse_complement(F[4])  
            gsm_key2 = F[0] + '\t' + ','.join(F[1:5])

            if gsm_key1 in key2gsm_info:
                gsm_info = key2gsm_info[gsm_key1]
                gsm_key = gsm_key1
            elif gsm_key2 in key2gsm_info: 
                gsm_info = key2gsm_info[gsm_key2] 
                gsm_key = gsm_key2
            else:
                continue
 
            for i in range(1, sig_num + 1):
                splice_pos_gsm_sig = gsm_info + '\t' + "TRUE" + '\t' + str(i)
                if splice_pos_gsm_sig not in splice_pos_gsm_sig2mem_sum: splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig] = 0
                splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig] = splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig] + float(F[header2ind["signature_" + str(i)]])
                if ctype == "STAD" and i == 1 and gsm_info == "donor" + '\t' + "5":
                    print '\t'.join(F)
                    print gsm_key + '\t' + F[header2ind["signature_" + str(i) ]]
  

    for splice_pos_gsm_sig in sorted(splice_pos_gsm_sig2mem_sum):
        splice_type, splice_pos, gsm, sig = splice_pos_gsm_sig.split('\t')
        cosm = sig2cosm[sig] if sig in sig2cosm else "NA"
        cor_value = sig2value[sig] if sig in sig2value else "NA"
        print >> hout_master, ctype + '\t' + splice_type + '\t' + splice_pos + '\t' + gsm + '\t' + str(sig) + '\t' + str(splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig]) + '\t' + cosm + '\t' + str(cor_value)

hout_master.close()


 

#! /usr/bin/env python

import sys, os, glob
import annot_utils.boundary
import subprocess

savnet_result_file = sys.argv[1]
output_file = sys.argv[2]
signature_dir = sys.argv[3]


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


annot_utils.boundary.make_boundary_info(output_file + ".refJunc.bed.gz", "hg19", False, "3,6", "6,1")

key2gsm = {}
with open(savnet_result_file, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = F[header2ind["Sample_Name"]] + '\t' + "chr" + F[header2ind["Mutation_Key"]]
        key2gsm[key] = 1
         


all_mem_files = glob.glob(signature_dir + "/*/mapping/mut_membership.txt")

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


    hout = open(output_file + ".tmp.bed", 'w')
    with open(mem_file, 'r') as hin:
        header = hin.readline()
        for line in hin:
            F = line.rstrip('\n').split('\t')
            print >> hout, F[1] + '\t' + str(int(F[2]) - 1) + '\t' + F[2]

    hout.close()

    hout = open(output_file + ".tmp.bed2", 'w')
    subprocess.call(["bedtools", "intersect", "-a", output_file + ".tmp.bed", "-b", output_file + ".refJunc.bed.gz", "-wb"], stdout = hout)
    hout.close()


    mut_key2splice_pos = {}
    with open(output_file + ".tmp.bed2", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mut_key = F[0] + '\t' + F[2]
            splice_pos = F[7] + '\t' + str(int(F[2]) - int(F[4])) if F[8] == "+" else F[7] + '\t' + str(int(F[5]) - int(F[2]) + 1)
            mut_key2splice_pos[mut_key] = splice_pos


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
            if F[1] + '\t' + F[2] not in mut_key2splice_pos: continue
            gsm_key1 = F[0] + '\t' + ','.join(F[1:5])
            F[3] = reverse_complement(F[3])
            F[4] = reverse_complement(F[4])  
            gsm_key2 = F[0] + '\t' + ','.join(F[1:5])
            is_gsm = "TRUE" if gsm_key1 in key2gsm or gsm_key2 in key2gsm else "FALSE"
 
            for i in range(1, sig_num + 1):
                splice_pos_gsm_sig = mut_key2splice_pos[F[1] + '\t' + F[2]] + '\t' + is_gsm + '\t' + str(i)
                if splice_pos_gsm_sig not in splice_pos_gsm_sig2mem_sum: splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig] = 0
                splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig] = splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig] + float(F[header2ind["signature_" + str(i)]])
                     

    for splice_pos_gsm_sig in sorted(splice_pos_gsm_sig2mem_sum):
        splice_type, splice_pos, gsm, sig = splice_pos_gsm_sig.split('\t')
        cosm = sig2cosm[sig] if sig in sig2cosm else "NA"
        cor_value = sig2value[sig] if sig in sig2value else "NA"
        print >> hout_master, ctype + '\t' + splice_type + '\t' + splice_pos + '\t' + gsm + '\t' + str(sig) + '\t' + str(splice_pos_gsm_sig2mem_sum[splice_pos_gsm_sig]) + '\t' + cosm + '\t' + str(cor_value)

hout_master.close()


subprocess.call(["rm", "-rf", output_file + ".refJunc.bed.gz"])
subprocess.call(["rm", "-rf", output_file + ".refJunc.bed.gz.tbi"])
subprocess.call(["rm", "-rf", output_file + ".tmp.bed"])
subprocess.call(["rm", "-rf", output_file + ".tmp.bed2"])

 

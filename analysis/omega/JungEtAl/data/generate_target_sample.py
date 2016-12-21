#! /usr/bin/env python

import glob

# all_files = glob.glob("../../output/*/mut_SJ_IR_list.txt")
all_files = glob.glob("../../gsm_out/gsm_file_list/*.mut_SJ_IR_list.txt")

sample_list_omega = {}

for list_file in sorted(all_files):
    with open(list_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "Sample_Name": continue
            sample_list_omega[F[0][:12]] = 1


hout = open("target_sample_list_JungEtAl.txt", 'w')
with open("sample_list_JungEtAl.txt", 'r') as hin:
    for line in hin:
        sample = line.rstrip('\n')
        # print sample[:12]
        if sample[:12] in sample_list_omega:
            print >> hout, sample

hout.close()



        
 

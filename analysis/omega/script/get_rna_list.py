#! /usr/bin/env python

import sys, os, glob

allfiles = glob.glob("../output_sv/*/sv_SJ_IR_Chimera_list.txt")

for gsm_input_file in sorted(allfiles):
    with open(gsm_input_file) as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "Sample_Name": continue
            print F[0]


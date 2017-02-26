#! /usr/bin/env python

import os, glob

all_intron_files = glob.glob("../intron_retention/*/TCGA-*.ir_simple_count.txt")


for intron_file in sorted(all_intron_files):
    sample = os.path.basename(intron_file).replace(".ir_simple_count.txt", "")
    if int(sample[13:15]) >= 10:
        print intron_file



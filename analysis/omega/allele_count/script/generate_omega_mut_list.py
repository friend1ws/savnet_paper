#! /usr/bin/env python

import sys, os, glob

output_file = sys.argv[1]
mut_dir = sys.argv[2]
bam_dir = sys.argv[3]

mut_files = glob.glob(mut_dir + "/*.genomon.genomon_mutation.result.filt.blacklist_filtered.txt")
bam_files = glob.glob(bam_dir + "/*/*.Aligned.sortedByCoord.out.bam")

# /home/omega3/omega_project/genomon/ACC/bam/TCGA-OR-A5KX-01/TCGA-OR-A5KX-01.markdup.bam


sample2mut_file = {}
for mut_file in sorted(mut_files):
    sample = os.path.basename(mut_file).replace(".genomon.genomon_mutation.result.filt.blacklist_filtered.txt", '')
    sample = sample[:15]
    if int(sample[13:15]) >= 10: continue
    sample2mut_file[sample] = mut_file

sample2bam_file = {}
for bam_file in sorted(bam_files):
    sample = os.path.basename(os.path.dirname(bam_file))
    sample = sample[:15]
    sample2bam_file[sample] = bam_file

hout = open(output_file, 'w')
for sample in sorted(sample2mut_file):
    if sample not in sample2bam_file: continue
    print >> hout, sample + '\t' + sample2mut_file[sample] + '\t' + sample2bam_file[sample] 

hout.close()




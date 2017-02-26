#! /usr/bin/env python

import sys, subprocess
import ebfilter.process_anno
import ebfilter.control_count

tumor_bam_file = sys.argv[1]
normal_bam_file = sys.argv[2]
anno_file = sys.argv[3]
output_file = sys.argv[4]
var_count_thres = 3
AF_tumor_thres = 0.05
AF_rel_thres = 10 

mapping_qual_thres = 20
base_qual_thres = 15
filter_flags = "UNMAP,SECONDARY,QCFAIL,DUP"
is_multi = False
is_loption = True 
region = ""

# generate pileup
ebfilter.process_anno.anno2pileup(anno_file, output_file + ".tmp.tumor.pileup", tumor_bam_file, 
                                  mapping_qual_thres, base_qual_thres, filter_flags, is_multi, is_loption, region)

ebfilter.process_anno.anno2pileup(anno_file, output_file + ".tmp.normal.pileup", normal_bam_file, 
                                  mapping_qual_thres, base_qual_thres, filter_flags, is_multi, is_loption, region)


pos2pileup_tumor = {}
with open(output_file + ".tmp.tumor.pileup", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        pos2pileup_tumor[F[0] + '\t' + F[1]] = [F[3], F[4], F[5]]

pos2pileup_normal = {}
with open(output_file + ".tmp.normal.pileup", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        pos2pileup_normal[F[0] + '\t' + F[1]] = [F[3], F[4], F[5]]


hout = open(output_file, 'w')

print >> hout, '\t'.join(["Chr", "Start", "End", "Ref", "Alt", "bases_tumor", "bases_normal", "misRate_tumor", "misRate_normal"])

with open(anno_file, 'r') as hin:
    for line in hin:
        chr, pos1, pos2, ref, alt, gene = line.rstrip('\n').split('\t')
        if alt == "-": pos1 = str(int(pos1) - 1)
        if chr + '\t' + pos1 not in pos2pileup_tumor: continue

        var = ""
        if ref != "-" and alt != "-":
            var = alt
        else:
            if ref == "-":
                var = "+" + alt
            elif alt == "-":
                var = "-" + ref

        if chr + '\t' + pos1 in pos2pileup_tumor:
            depth, base_bar, qual_bar = pos2pileup_tumor[chr + '\t' + pos1]
            var_tumor_p, depth_tumor_p, var_tumor_n, depth_tumor_n = ebfilter.control_count.varCountCheck(var, depth, base_bar, qual_bar, base_qual_thres, False)
        else:
            var_tumor_p, depth_tumor_p, var_tumor_n, depth_tumor_n = 0, 0, 0, 0

        if chr + '\t' + pos1 in pos2pileup_normal:
            depth, base_bar, qual_bar = pos2pileup_normal[chr + '\t' + pos1]
            var_normal_p, depth_normal_p, var_normal_n, depth_normal_n = ebfilter.control_count.varCountCheck(var, depth, base_bar, qual_bar, base_qual_thres, False)
        else:
            var_normal_p, depth_normal_p, var_normal_n, depth_normal_n = 0, 0, 0, 0

        if depth_tumor_p + depth_tumor_n == 0: continue
        if depth_normal_p + depth_normal_n == 0: continue

        AF_tumor = float(var_tumor_p + var_tumor_n) / float(depth_tumor_p + depth_tumor_n)
        AF_normal = float(var_normal_p + var_normal_n) / float(depth_normal_p + depth_normal_n)

        if var_tumor_p + var_tumor_n < var_count_thres: continue
        if AF_tumor < AF_tumor_thres: continue

        if AF_normal > 0 and AF_tumor / AF_normal < AF_rel_thres: continue

        print >> hout, '\t'.join([chr, pos1, pos2, ref, alt]) + '\t' + ','.join([str(x) for x in [depth_tumor_p, var_tumor_p, depth_tumor_n, var_tumor_n]]) + '\t' + \
              ','.join([str(x) for x in [depth_normal_p, var_normal_p, depth_normal_n, var_normal_n]]) + '\t' + str(AF_tumor) + '\t' + str(AF_normal)
           
hout.close()


subprocess.call(["rm", "-rf", output_file + ".tmp.tumor.pileup"])
subprocess.call(["rm", "-rf", output_file + ".tmp.normal.pileup"])



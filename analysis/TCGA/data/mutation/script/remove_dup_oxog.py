#! /usr/bin/env python

import sys, subprocess
from operator import itemgetter
import pysam

mut_file = sys.argv[1]
oxog_file = sys.argv[2]
canonical_file = sys.argv[3]
noncanonical_file = sys.argv[4]

canonical_db = pysam.TabixFile(canonical_file)
noncanonical_db = pysam.TabixFile(noncanonical_file)


key2oxog = {}
with open(oxog_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')

        oxog_flag = 0
        if F[header2ind["Ref"]] == "C" and F[header2ind["Alt"]] == "A":
            if int(F[header2ind["alt_F1R2"]]) < 2: oxog_flag = 1
            if float(F[header2ind["alt_F2R1"]]) / (float(F[header2ind["alt_F1R2"]]) + float(F[header2ind["alt_F2R1"]])) >= 0.9: oxog_flag = 1
        elif F[header2ind["Ref"]] == "G" and F[header2ind["Alt"]] == "T":
            if int(F[header2ind["alt_F2R1"]]) < 2: oxog_flag = 1
            if float(F[header2ind["alt_F1R2"]]) / (float(F[header2ind["alt_F1R2"]]) + float(F[header2ind["alt_F2R1"]])) >= 0.9: oxog_flag = 1

        if oxog_flag == 1:
            key = '\t'.join([ F[header2ind[x]] for x in ["Chr", "Start", "End", "Ref", "Alt"] ])
            key2oxog[key] = 1



def check_junc(mut):

    # check overlap with cannical splicing motif sites
    tabixErrorFlag = 0
    try:
        record = canonical_db.fetch(mut[0], int(mut[1]) - 10, int(mut[2]) + 10)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0 and record is not None:
        for record_line in record:
            junction = record_line.split('\t')
            if mut[2] >= int(junction[1]) + 1 and mut[1] <= int(junction[2]):
                return 2

    # check overlap with non-canonical splicing motif sites
    tabixErrorFlag = 0
    try:
        record = noncanonical_db.fetch(mut[0], int(mut[1]) - 10, int(mut[2]) + 10)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    junc = 0
    if tabixErrorFlag == 0 and record is not None:
        for record_line in record:
            junction = record_line.split('\t')
            if mut[2] >= int(junction[1]) + 1 and mut[1] <= int(junction[2]):
                return 1

    # no overlap
    return 0


def check_indel(mut):
    if len(mut[3]) > 1: # deletion
        return 2
    elif len(mut[4]) > 1: # insertion
        return 1 # SNV
    else:
        return 0

def mut_compare(mut1, mut2):
    junc1 = check_junc(mut1)
    junc2 = check_junc(mut2)

    print >> sys.stderr, '\t'.join([str(x) for x in mut1]) + '\t' + str(junc1)
    print >> sys.stderr, '\t'.join([str(x) for x in mut2]) + '\t' + str(junc2)

    if junc1 > junc2:
        return mut1
    elif junc2 > junc1:
        return mut2
    else:

        is_indel1 = check_indel(mut1)
        is_indel2 = check_indel(mut2)

        if is_indel1 > is_indel2:
            return mut1
        elif is_indel2 > is_indel1:
            return mut2
        else:

            if mut1[1] >= mut2[1]:
                return mut1
            else:
                return mut2



mut_list = []
header2ind = {}
with open(mut_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0].startswith("#"): continue 
        if F[0] == "Chr": 
            for (i, cname) in enumerate(F):
                header2ind[cname] = i

            print '\t'.join(F[:10])
            continue

        key = '\t'.join([F[header2ind[x]] for x in ["Chr", "Start", "End", "Ref", "Alt"]])
        if key in key2oxog: continue

        if float(F[header2ind["misRate_tumor"]]) < 0.05: continue
        mut_list.append((F[0], int(F[1]), int(F[2]), F[3], F[4], F[5], F[6], F[7], F[8], F[9]))



temp_chr = ""
temp_start = 0
temp_end = 0
temp_mut = []
for mut in sorted(mut_list, key = itemgetter(0, 1, 2)):

    if str(mut[1]) == "116412035":
        pass

    if mut[0] != temp_chr or mut[1] > temp_end + 5:
        if temp_chr != "":
            print '\t'.join([str(x) for x in temp_mut])
        temp_mut = mut
    else:
        # check duplicate
        temp_mut = mut_compare(temp_mut, mut)

    temp_chr = temp_mut[0]
    temp_start = temp_mut[1]
    temp_end = temp_mut[2]

if temp_chr != "":
    print '\t'.join([str(x) for x in temp_mut])


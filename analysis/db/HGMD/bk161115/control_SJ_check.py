#! /usr/bin/env python

import sys, gzip
import pysam

hgmd_cs_db = pysam.TabixFile("HGMD.CS.bed.gz")

with gzip.open("/home/yshira/project/inframe_junc/output/control.bed.gz", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        if F[0] in ["MT", "Y", "hs37d5"]: continue
        if F[0].startswith("GL0"): continue

        tabixErrorFlag = 0
        try:
            records = hgmd_cs_db.fetch("chr" + F[0], int(F[1]) - 30, int(F[2]) + 30)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

        mut_keys = []
        genes = []
        for record_line in records:
            record = record_line.split('\t')
            record_info = record[5].split(';')
            mut_keys.append(record[0] + ':' + record[1] + '-' + record[2])
            genes.append(record_info[4])

        mut_keys = list(set(mut_keys))
        genes = list(set(genes))

        if len(mut_keys) > 0:
            print '\t'.join(F) + '\t' + ','.join(genes) + '\t' + ','.join(mut_keys)




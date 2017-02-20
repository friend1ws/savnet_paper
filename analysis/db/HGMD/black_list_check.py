#! /usr/bin/env python

import sys, gzip
import pysam


black_list = "/home/kchiba/work_directory/work_hotspot/map_reduce/omega_project.black_list.bed.gz"
hgmd = "HGMD.CS.bed.gz"

black_db = pysam.TabixFile(black_list)


with gzip.open(hgmd) as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        tchr = F[0].replace('chr', "")
        FF = F[5].split(';')

        tabixErrorFlag = 0
        try:
            records = black_db.fetch(tchr, int(F[2]) - 10, int(F[2]) + 10)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

        mut_keys = []
        genes = []
        for record_line in records:
            record = record_line.split('\t')
            if F[1] == record[1] and F[2] == record[2] and F[3] == record[3] and F[4] == record[4]:
                print '\t'.join(F) + '\t' + FF[4] + '\t' + FF[12] + '\t' + record_line
 

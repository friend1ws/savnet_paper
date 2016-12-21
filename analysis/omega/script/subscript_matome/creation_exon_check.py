#! /usr/bin/env python

import sys, os, pysam, annot_utils.exon

input_file = sys.argv[1]
output_file = sys.argv[2]

annot_utils.exon.make_exon_info(output_file + ".exon.bed.gz", "ref", "hg19", True, False)
exon_info_tb = pysam.TabixFile(output_file + ".exon.bed.gz")

hout = open(output_file, 'w')
header2ind = {}
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    print >> hout, '\t'.join(header) + '\t' + "Exon_Strand" + '\t' + "Mut_Pos" + '\t' + "Exon_Start" + '\t' + "Exon_End"
 
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[header2ind["Motif_Type"]] in ["splicing donor disruption", "splicing acceptor disruption"]: continue

        FF = F[header2ind["Mutation_Key"]].split(',')

        tabixErrorFlag = 0
        try:
            records = exon_info_tb.fetch(FF[0], int(FF[1]) - 5, int(FF[1]) + 5)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

        exon_start = []
        exon_end = []
        exon_strand = []
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                exon_start.append(record[1])
                exon_end.append(record[2])
                exon_strand.append(record[5])

        if len(exon_start) == 0: continue

        exon_start = list(set(exon_start))
        exon_end = list(set(exon_end))
        exon_strand = list(set(exon_strand))

        if len(exon_start) > 1 or len(exon_end) > 1 or len(exon_strand) > 1: continue

        print >> hout, '\t'.join(F) + '\t' + exon_strand[0] + '\t' + FF[1] + '\t' + exon_start[0] + '\t' + exon_end[0]


hout.close()

os.unlink(output_file + ".exon.bed.gz")
os.unlink(output_file + ".exon.bed.gz.tbi")


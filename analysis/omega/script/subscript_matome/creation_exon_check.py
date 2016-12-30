#! /usr/bin/env python

import sys, os, pysam, re, annot_utils.exon

input_file = sys.argv[1]
output_file = sys.argv[2]

annot_utils.exon.make_exon_info(output_file + ".exon.bed.gz", "refseq", "hg19", True, False)
exon_info_tb = pysam.TabixFile(output_file + ".exon.bed.gz")

hout = open(output_file, 'w')
header2ind = {}
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    print >> hout, '\t'.join(header) + '\t' + "Exon_Strand" + '\t' + "Junc_Pos" + '\t' + "Exon_Start" + '\t' + "Exon_End"
 
    for line in hin:
        F = line.rstrip('\n').split('\t')
        # if F[header2ind["Mutation_Type"]] in ["splicing donor disruption", "splicing acceptor disruption"]: continue
        if F[header2ind["Splicing_Class"]] not in ["alternative-3'-splice-site", "alternative-5'-splice-site", 
                                                 "intronic-alternative-3'-splice-site", "intronic-alternative-5'-splice-site"]: continue

        if F[header2ind["Mutation_Type"]] == "splicing donor disruption" and F[header2ind["Splicing_Class"]].endswith("3'-splice-site"): continue
        if F[header2ind["Mutation_Type"]] == "splicing acceptor disruption" and F[header2ind["Splicing_Class"]].endswith("5'-splice-site"): continue

        FF = F[header2ind["Mutation_Key"]].split(',')

        junc_match = re.match(r'([^ \t\n\r\f\v,]+)\:(\d+)\-(\d+)', F[header2ind["Splicing_Key"]])
        if junc_match is None:
            print "Something is wrong"
            sys.exit(1)
        junc_chr, junc_start, junc_end = junc_match.group(1), junc_match.group(2), junc_match.group(3)

        tabixErrorFlag = 0
        try:
            records = exon_info_tb.fetch(FF[0], int(FF[1]) - 30, int(FF[1]) + 30)
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

        if F[header2ind["Mutation_Key"]] == "1,27100390,G,A":
            pass

        if exon_strand[0] == "+" and F[header2ind["Splicing_Class"]].endswith("5'-splice-site") or \
           exon_strand[0] == "-" and F[header2ind["Splicing_Class"]].endswith("3'-splice-site"):
            print >> hout, '\t'.join(F) + '\t' + exon_strand[0] + '\t' + junc_start + '\t' + exon_start[0] + '\t' + exon_end[0]
        else:
            print >> hout, '\t'.join(F) + '\t' + exon_strand[0] + '\t' + junc_end + '\t' + exon_start[0] + '\t' + exon_end[0]

hout.close()

os.unlink(output_file + ".exon.bed.gz")
os.unlink(output_file + ".exon.bed.gz.tbi")


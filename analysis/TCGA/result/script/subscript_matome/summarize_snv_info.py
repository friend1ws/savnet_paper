#! /usr/bin/env python

import sys, re
# from subscript_matome import my_seq
import my_utils.seq
import pysam

input_file = sys.argv[1]
output_file = sys.argv[2]
reference = sys.argv[3]
# hgmd_file = sys.argv[4]
spidex_file = sys.argv[4]


key2exists = {}
header2ind = {}

seq_margin = 100

# hgmd_db = pysam.TabixFile(hgmd_file)
spidex_db = pysam.TabixFile(spidex_file)

hout = open(output_file, 'w')
print >> hout, '\t'.join(["Cancer_Type", "Sample_Name", "Gene_Symbol", "Mutation_Key", "Motif_Pos", 
                          "Motif_Seq", "Rel_Pos", "Ref_Base", "Alt_Base", "Mutation_Type", "Is_Canonical", "SPIDEX"])
 
with open(input_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    for i in range(len(header)):
        header2ind[header[i]] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        # if F[header2ind["Mutation_Type"]] not in ["splicing donor disruption", "splicing acceptor disruption"]: continue
        key = F[header2ind["Cancer_Type"]] + '\t' + F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Gene_Symbol"]] + '\t' + \
              F[header2ind["Mutation_Key"]] + '\t' + F[header2ind["Motif_Pos"]]

        if key in key2exists: continue


        mchr, mpos, mref, malt = F[header2ind["Mutation_Key"]].split(',')
        if len(mref) > 1 or len(malt) > 1: continue

        """
        ##########
        # check hgmd
        tabixErrorFlag = 0
        try:
            records = hgmd_db.fetch("chr" + mchr, int(mpos) - 3, int(mpos) + 3)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

        hgmd_info = []
        for record_line in records:
            record = record_line.split('\t')
            record_info = record[5].split(';')
            if mpos == record[2]:
                hgmd_info.append(record_info[-1].strip('"'))

        hgmd_info_print = ';'.join(list(set(hgmd_info))) if len(hgmd_info) > 0 else "---"
        ##########
        """

        ##########
        # check spidex 
        tabixErrorFlag = 0
        try:
            records = spidex_db.fetch(mchr, int(mpos) - 3, int(mpos) + 3)
        except Exception as inst:
            print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1

        spidex = "---" 
        for record_line in records:
            record = record_line.split('\t')
            record_info = record[5].split(';')
            if mpos == record[2] and mref == record[3] and malt == record[4]:
                spidex = record[6]

        motif_pos, strand = F[header2ind["Motif_Pos"]].split(',')

        pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', motif_pos)
        schr, sstart, send = pos_match.group(1), pos_match.group(2), pos_match.group(3)

        if strand == "+":
            motif_seq = my_utils.seq.get_seq(reference, motif_pos)
            rel_pos = int(mpos) - int(sstart) + 1
            print >> hout, key + '\t' + motif_seq + '\t' + str(rel_pos) + '\t' + mref + '\t' + malt + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]] + '\t' + spidex 

        else:
            motif_seq = my_utils.seq.reverse_complement(my_utils.seq.get_seq(reference, motif_pos))
            rel_pos = int(send) - int(mpos) + 1
            print >> hout, key + '\t' + motif_seq + '\t' + str(rel_pos) + '\t' + my_utils.seq.reverse_complement(mref) + '\t' + my_utils.seq.reverse_complement(malt) + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]] + '\t' + spidex 

        key2exists[key] = 1

hout.close()



#! /usr/bin/env python

import sys, re
# from subscript_matome import my_seq
import my_seq

input_file = sys.argv[1]
output_file = sys.argv[2]
reference = sys.argv[3]


key2exists = {}
header2ind = {}

seq_margin = 100

hout = open(output_file, 'w')
print >> hout, '\t'.join(["Cancer_Type", "Sample_Name", "Gene_Symbol", "Mutation_Key", 
                          "Motif_Seq", "Rel_Pos", "Ref_Base", "Alt_Base", "Motif_Type", "Is_Canonocal"])
 
with open(input_file, 'r') as hin:

    header = hin.readline().rstrip('\n').split('\t')
    for i in range(len(header)):
        header2ind[header[i]] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        # if F[header2ind["Mutation_Type"]] not in ["splicing donor disruption", "splicing acceptor disruption"]: continue
        key = F[header2ind["Cancer_Type"]] + '\t' + F[header2ind["Sample_Name"]] + '\t' + F[header2ind["Gene_Symbol"]] + '\t' + F[header2ind["Mutation_Key"]]
        if key in key2exists: continue


        mchr, mpos, mref, malt = F[header2ind["Mutation_Key"]].split(',')
        if len(mref) > 1 or len(malt) > 1: continue

        motif_pos, strand = F[header2ind["Motif_Pos"]].split(',')

        pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', motif_pos)
        schr, sstart, send = pos_match.group(1), pos_match.group(2), pos_match.group(3)

        if strand == "+":
            motif_seq = my_seq.get_seq(reference, motif_pos)
            rel_pos = int(mpos) - int(sstart) + 1
            print >> hout, key + '\t' + motif_seq + '\t' + str(rel_pos) + '\t' + mref + '\t' + malt + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]]

        else:
            motif_seq = my_seq.reverse_complement(my_seq.get_seq(reference, motif_pos))
            rel_pos = int(send) - int(mpos) + 1
            print >> hout, key + '\t' + motif_seq + '\t' + str(rel_pos) + '\t' + my_seq.reverse_complement(mref) + '\t' + my_seq.reverse_complement(malt) + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]]

        key2exists[key] = 1

hout.close()



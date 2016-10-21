#! /usr/bin/env python

import sys, re
# from subscript_matome import my_seq
import my_seq
import pysam
from annot_utils.junction import *

input_file = sys.argv[1]
output_file = sys.argv[2]
reference = sys.argv[3]

# generate junction bed file
# make_junction_info(output_file + ".refJunc.bed.gz", "hg19", True, "1,0", "0,1")
ref_junc_tb = pysam.TabixFile(output_file + ".refJunc.bed.gz")

key2exists = {}
header2ind = {}

seq_margin = 100

hout = open(output_file, 'w')
print >> hout, '\t'.join(["Cancer_Type", "Sample_Name", "Gene_Symbol", "Mutation_Key", "Motif_Pos", 
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


        mut_chr, mut_pos, mut_ref, mut_alt = F[header2ind["Mutation_Key"]].split(',')
        if len(mut_ref) > 1 or len(mut_alt) > 1: continue

        pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', F[header2ind["Splicing_Key"]])
        sp_chr, sp_start, sp_end = pos_match.group(1), pos_match.group(2), pos_match.group(3)

        motif_pos, motif_strand = F[header2ind["Motif_Pos"]].split(',')

        pos_match = re.match(r'([\w\d]+)\:(\d+)\-(\d+)', motif_pos)
        motif_chr, motif_start, motif_end = pos_match.group(1), pos_match.group(2), pos_match.group(3)

        motif_bias = 0
        ##########
        # for splicing motif creating mutations, modify the motif position for alignment redundancy
        if F[header2ind["Mutation_Type"]] in ["splicing donor creation", "splicing acceptor creation"]:

            if motif_strand == "+" and F[header2ind["Mutation_Type"]] == "splicing donor creation" or \
                motif_strand == "-" and F[header2ind["Mutation_Type"]] == "splicing acceptor creation":
                new_junc, annot_junc = sp_start, sp_end
            else:
                new_junc, annot_junc = sp_end, sp_start

            ####################
            # get the records for control junction data for the current position
            tabixErrorFlag = 0
            try:
                records = ref_junc_tb.fetch(sp_chr, int(annot_junc) - 5, int(annot_junc) + 5)
            except Exception as inst:
                print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorMsg = str(inst.args)
                tabixErrorFlag = 1
    
            ####################
            # for each record in control junction extracted, check the consistency with the current junction
            if tabixErrorFlag == 0:
               
                for record_line in records:
                    record = record_line.split('\t')
   
                    if record[3] != F[header2ind["Gene_Symbol"]]: continue 
                    if record[4] == "donor" and F[header2ind["Mutation_Type"]] == "splicing donor creation": continue
                    if record[4] == "acceptor" and F[header2ind["Mutation_Type"]] == "splicing acceptor creation": continue
 
                    if record[4] == "donor" and record[5] == "+" or record[4] == "acceptor" and record[5] == "-":  
                        annot_junc_true = int(record[2]) + 1
                        # motif_bias = int(annot_junc_true) - int(annot_junc)
                    else:
                        annot_junc_true = int(record[2]) - 1
                        # motif_bias = int(annot_junc) - int(annot_junc_true)

                    motif_bias = int(annot_junc_true) - int(annot_junc)
                    """
                    if int(annot_junc) != int(annot_junc_true): 
                        # print str(annot_junc) + '\t' + str(annot_junc_true)
                        # print record_line
                        # print '\t'.join(F)

                        if motif_strand == "+" and F[header2ind["Mutation_Type"]] == "splicing donor creation" or \
                            motif_strand == "-" and F[header2ind["Mutation_Type"]] == "splicing acceptor creation":
                                motif_bias = int(annot_junc_true) - int(annot_junc)
                        else:
                                motif_bias = int(annot_junc) - int(annot_junc_true)
                    """


        if motif_bias != 0:
            motif_start, motif_end = str(int(motif_start) + motif_bias), str(int(motif_end) + motif_bias)
            motif_pos = motif_chr + ':' +  motif_start + '-' + motif_end           

        # if motif_bias == 0: continue

        if F[header2ind["Mutation_Key"]] in ["14,89077251,G,A", "15,42127222,G,A", "1,155932910,C,T"]:
            print motif_bias
            print motif_pos
            print F[header2ind["Motif_Pos"]]
        ##########
 
        if motif_strand == "+":
            motif_seq = my_seq.get_seq(reference, motif_pos)
            rel_pos = int(mut_pos) - int(motif_start) + 1
            print >> hout, key + '\t' + motif_pos + '\t' + motif_seq + '\t' + str(rel_pos) + '\t' + mut_ref + '\t' + mut_alt + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]]

        else:
            motif_seq = my_seq.reverse_complement(my_seq.get_seq(reference, motif_pos))
            rel_pos = int(motif_end) - int(mut_pos) + 1
            print >> hout, key + '\t' + motif_pos + '\t' + motif_seq + '\t' + str(rel_pos) + '\t' + my_seq.reverse_complement(mut_ref) + '\t' + my_seq.reverse_complement(mut_alt) + '\t' + F[header2ind["Mutation_Type"]] + '\t' + F[header2ind["Is_Canonical"]]

        key2exists[key] = 1

hout.close()


